#!/bin/bash
# Check getfastq workflow status and progress with enhanced tracking
#
# Usage:
#   bash scripts/rna/check_getfastq_status.sh [--species SPECIES] [--verbose] [--help]
#
# Options:
#   --species SPECIES    Species name (default: pogonomyrmex_barbatus)
#   --verbose, -v        Show detailed sample lists for each category
#   --help, -h            Show this help message
#
# Examples:
#   bash scripts/rna/check_getfastq_status.sh
#   bash scripts/rna/check_getfastq_status.sh --species camponotus_floridanus
#   bash scripts/rna/check_getfastq_status.sh --verbose

set -e

# Parse command line arguments
VERBOSE=0
SPECIES="pogonomyrmex_barbatus"

while [[ $# -gt 0 ]]; do
    case $1 in
        --species)
            SPECIES="$2"
            shift 2
            ;;
        --verbose|-v)
            VERBOSE=1
            shift
            ;;
        --help|-h)
            echo "Check getfastq workflow status and progress"
            echo ""
            echo "Usage:"
            echo "  $0 [--species SPECIES] [--verbose] [--help]"
            echo ""
            echo "Options:"
            echo "  --species SPECIES    Species name (default: pogonomyrmex_barbatus)"
            echo "  --verbose, -v        Show detailed sample lists for each category"
            echo "  --help, -h            Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0"
            echo "  $0 --species camponotus_floridanus"
            echo "  $0 --verbose"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Use --help for usage information" >&2
            exit 1
            ;;
    esac
done

FASTQ_DIR="output/amalgkit/${SPECIES}/fastq/getfastq"
QUANT_DIR="output/amalgkit/${SPECIES}/quant"
LOG_DIR="output/amalgkit/${SPECIES}/logs"
METADATA="output/amalgkit/${SPECIES}/work/metadata/metadata.tsv"

# Color codes (with fallback to no color if not supported)
if [ -t 1 ] && command -v tput >/dev/null 2>&1; then
    RED=$(tput setaf 1 2>/dev/null || echo "")
    GREEN=$(tput setaf 2 2>/dev/null || echo "")
    YELLOW=$(tput setaf 3 2>/dev/null || echo "")
    BLUE=$(tput setaf 4 2>/dev/null || echo "")
    RESET=$(tput sgr0 2>/dev/null || echo "")
else
    RED=""
    GREEN=""
    YELLOW=""
    BLUE=""
    RESET=""
fi

echo "═══════════════════════════════════════════════════════════════════════════════"
echo "📊 RNA-seq Workflow Status Check"
echo "═══════════════════════════════════════════════════════════════════════════════"
echo "Species: ${BLUE}$SPECIES${RESET}"
echo ""

# Check if directories exist
if [ ! -d "$FASTQ_DIR" ]; then
    echo "${YELLOW}⚠️  FASTQ directory not found: $FASTQ_DIR${RESET}"
    echo "   This may indicate the workflow hasn't started yet, or the species name is incorrect."
    echo "   Use --species to specify a different species, or check the directory structure."
    echo ""
    # Don't exit - continue to show what we can
fi

# Get total samples from metadata
if [ -f "$METADATA" ]; then
    TOTAL_SAMPLES=$(tail -n +2 "$METADATA" | wc -l)
    echo "📋 Total samples in metadata: $TOTAL_SAMPLES"
else
    echo "${YELLOW}⚠️  Metadata file not found: $METADATA${RESET}"
    TOTAL_SAMPLES=0
fi

echo ""

# Get active processes and extract sample IDs
declare -A ACTIVE_PROCESS_SAMPLES
declare -A PROCESS_PIDS
ACTIVE_QUANT_SAMPLES=()

if command -v ps >/dev/null 2>&1; then
    # Get getfastq processes with sample IDs
    while IFS= read -r line; do
        pid=$(echo "$line" | awk '{print $2}')
        # Extract sample ID from metadata path
        sample=$(echo "$line" | grep -o "metadata\.download\.SRR[0-9]*" | sed 's/metadata\.download\.//' | head -1)
        if [ -n "$sample" ] && [ -n "$pid" ]; then
            ACTIVE_PROCESS_SAMPLES["$sample"]="$pid"
            PROCESS_PIDS["$pid"]="$sample"
        fi
    done < <(ps aux | grep -E "amalgkit getfastq" | grep -v grep | awk '{pid=$2; cmd=""; for(i=11;i<=NF;i++) cmd=cmd" "$i; print pid, cmd}' | sort -u -k1,1 || true)
    
    # Get quant processes with sample IDs
    while IFS= read -r line; do
        # Extract sample ID from quant log file names in process args
        sample=$(echo "$line" | grep -o "SRR[0-9]*" | head -1)
        if [ -n "$sample" ]; then
            ACTIVE_QUANT_SAMPLES+=("$sample")
        fi
    done < <(ps aux | grep -E "amalgkit quant" | grep -v grep || true)
fi

# Count sample directories
TOTAL_DIRS=$(find "$FASTQ_DIR" -mindepth 1 -maxdepth 1 -type d -name "SRR*" 2>/dev/null | wc -l)

# Count by status and track samples
CONVERTED_NOT_QUANTIFIED=0
CONVERTED_QUANTIFIED=0
CONVERTED_QUANTIFIED_DELETED=0
QUANTIFYING=0
CONVERTING=0
DOWNLOADING=0

# Track samples for each category
declare -A SAMPLE_STATUS
declare -A SAMPLE_CATEGORIES
declare -A ACTIVE_DOWNLOADS
declare -A STUCK_SAMPLES
declare -A SAMPLE_LAST_ACTIVITY

# Initialize category arrays
COMPLETED_SAMPLES=()
QUANTIFIED_SAMPLES=()
QUANTIFYING_SAMPLES=()
NOT_QUANTIFIED_SAMPLES=()
CONVERTING_SAMPLES=()
DOWNLOADING_SAMPLES=()

CURRENT_TIME=$(date +%s)

# Check log files for active downloads and track timestamps
if [ -d "$LOG_DIR" ]; then
    for log_file in "$LOG_DIR"/*getfastq*.stdout.log "$LOG_DIR"/*getfastq*.log; do
        if [ -f "$log_file" ]; then
            # Check if log was modified in last 5 minutes
            MOD_TIME=$(stat -c %Y "$log_file" 2>/dev/null || stat -f %m "$log_file" 2>/dev/null || echo "0")
            if [ "$MOD_TIME" != "0" ]; then
                TIME_DIFF=$((CURRENT_TIME - MOD_TIME))
                if [ "$TIME_DIFF" -lt 300 ]; then
                    # Extract sample ID from log filename
                    sample=$(basename "$log_file" | grep -o "SRR[0-9]*" | head -1)
                    if [ -n "$sample" ]; then
                        ACTIVE_DOWNLOADS["$sample"]=1
                        SAMPLE_LAST_ACTIVITY["$sample"]=$MOD_TIME
                    fi
                fi
            fi
        fi
    done
fi

# Check each sample directory
for sample_dir in "$FASTQ_DIR"/SRR*/; do
    if [ -d "$sample_dir" ]; then
        sample=$(basename "$sample_dir")
        has_fastq=$(find "$sample_dir" -name "*.fastq*" ! -name "*.sra" 2>/dev/null | wc -l)
        has_sra=$(find "$sample_dir" -name "*.sra" 2>/dev/null | wc -l)
        is_quantified=$([ -f "$QUANT_DIR/$sample/abundance.tsv" ] && echo "1" || echo "0")
        is_quantifying=0
        
        # Get directory modification time for stuck detection
        DIR_MOD_TIME=$(stat -c %Y "$sample_dir" 2>/dev/null || stat -f %m "$sample_dir" 2>/dev/null || echo "0")
        if [ "$DIR_MOD_TIME" != "0" ]; then
            SAMPLE_LAST_ACTIVITY["$sample"]=$DIR_MOD_TIME
        fi
        
        # Check if this sample is currently being quantified
        for quant_sample in "${ACTIVE_QUANT_SAMPLES[@]}"; do
            if [ "$quant_sample" = "$sample" ]; then
                is_quantifying=1
                break
            fi
        done
        
        # Check for stuck samples
        if [ "$DIR_MOD_TIME" != "0" ]; then
            TIME_SINCE_ACTIVITY=$((CURRENT_TIME - DIR_MOD_TIME))
            # Stuck if downloading > 30 min or converting > 1 hour
            if [ "$has_sra" -gt 0 ] && [ "$has_fastq" -eq 0 ] && [ "$TIME_SINCE_ACTIVITY" -gt 3600 ]; then
                STUCK_SAMPLES["$sample"]="converting_stuck"
            elif [ "$has_sra" -eq 0 ] && [ "$has_fastq" -eq 0 ] && [ "$TIME_SINCE_ACTIVITY" -gt 1800 ]; then
                STUCK_SAMPLES["$sample"]="downloading_stuck"
            fi
        fi
        
        if [ "$is_quantifying" -eq 1 ]; then
            QUANTIFYING=$((QUANTIFYING + 1))
            SAMPLE_STATUS["$sample"]="quantifying"
            QUANTIFYING_SAMPLES+=("$sample")
        elif [ "$has_fastq" -gt 0 ]; then
            if [ "$is_quantified" -eq 1 ]; then
                CONVERTED_QUANTIFIED=$((CONVERTED_QUANTIFIED + 1))
                SAMPLE_STATUS["$sample"]="converted_quantified"
                QUANTIFIED_SAMPLES+=("$sample")
            else
                CONVERTED_NOT_QUANTIFIED=$((CONVERTED_NOT_QUANTIFIED + 1))
                SAMPLE_STATUS["$sample"]="converted_not_quantified"
                NOT_QUANTIFIED_SAMPLES+=("$sample")
            fi
        elif [ "$has_sra" -gt 0 ]; then
            CONVERTING=$((CONVERTING + 1))
            SAMPLE_STATUS["$sample"]="converting"
            CONVERTING_SAMPLES+=("$sample")
        else
            # Check if actively downloading based on log activity or process
            if [ -n "${ACTIVE_DOWNLOADS[$sample]}" ] || [ -n "${ACTIVE_PROCESS_SAMPLES[$sample]}" ]; then
                DOWNLOADING=$((DOWNLOADING + 1))
                SAMPLE_STATUS["$sample"]="downloading"
                DOWNLOADING_SAMPLES+=("$sample")
            else
                # Empty directory but no recent activity - might be stuck or not started
                DOWNLOADING=$((DOWNLOADING + 1))
                SAMPLE_STATUS["$sample"]="downloading"
                DOWNLOADING_SAMPLES+=("$sample")
            fi
        fi
    fi
done

# Check for quantified samples without FASTQ directories (deleted)
if [ -d "$QUANT_DIR" ]; then
    for quant_sample_dir in "$QUANT_DIR"/SRR*/; do
        if [ -d "$quant_sample_dir" ]; then
            sample=$(basename "$quant_sample_dir")
            if [ -f "$quant_sample_dir/abundance.tsv" ]; then
                # Check if FASTQ directory exists or has FASTQ files
                if [ ! -d "$FASTQ_DIR/$sample" ] || [ $(find "$FASTQ_DIR/$sample" -name "*.fastq*" ! -name "*.sra" 2>/dev/null | wc -l) -eq 0 ]; then
                    # Only count if not already counted
                    if [ -z "${SAMPLE_STATUS[$sample]}" ]; then
                        CONVERTED_QUANTIFIED_DELETED=$((CONVERTED_QUANTIFIED_DELETED + 1))
                        SAMPLE_STATUS["$sample"]="completed"
                        COMPLETED_SAMPLES+=("$sample")
                    fi
                fi
            fi
        fi
    done
fi

# Calculate not started: total - all categories that have been started
STARTED=$((CONVERTED_NOT_QUANTIFIED + CONVERTED_QUANTIFIED + CONVERTED_QUANTIFIED_DELETED + QUANTIFYING + CONVERTING + DOWNLOADING))
NOT_STARTED=$((TOTAL_SAMPLES - STARTED))

# Progress bar function
draw_progress_bar() {
    local current=$1
    local total=$2
    local width=40
    local percent=0
    
    if [ "$total" -gt 0 ]; then
        percent=$((current * 100 / total))
    fi
    
    local filled=$((current * width / total))
    local empty=$((width - filled))
    
    printf "["
    if [ "$filled" -gt 0 ]; then
        printf "%${filled}s" | tr ' ' '='
    fi
    if [ "$empty" -gt 0 ]; then
        printf "%${empty}s" | tr ' ' '-'
    fi
    printf "] %3d%%" "$percent"
}

# Summary section at top
echo "📊 Quick Summary:"
if [ "$TOTAL_SAMPLES" -gt 0 ]; then
    PERCENT_COMPLETE=$((CONVERTED_QUANTIFIED_DELETED * 100 / TOTAL_SAMPLES))
    PERCENT_STARTED=$((STARTED * 100 / TOTAL_SAMPLES))
    echo "   ${GREEN}✅ Completed (goal state):${RESET}     $CONVERTED_QUANTIFIED_DELETED / $TOTAL_SAMPLES ($PERCENT_COMPLETE%)"
    echo "   ${BLUE}🔄 In Progress:${RESET}                $((STARTED - CONVERTED_QUANTIFIED_DELETED)) / $TOTAL_SAMPLES"
    echo "   ${YELLOW}⏸️  Not Started:${RESET}                $NOT_STARTED / $TOTAL_SAMPLES"
    echo ""
    echo "   Progress: $(draw_progress_bar $CONVERTED_QUANTIFIED_DELETED $TOTAL_SAMPLES)"
fi
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""

echo "📁 Detailed Status Breakdown:"
echo "   ${GREEN}✅ Converted, quantified, deleted:${RESET} $CONVERTED_QUANTIFIED_DELETED (goal state - RNA quantified, FASTQ deleted)"
echo "   ${BLUE}🔬 Converted and quantified:${RESET}      $CONVERTED_QUANTIFIED (FASTQ still present, can be deleted)"
echo "   ${BLUE}🧬 Quantifying:${RESET}                   $QUANTIFYING (currently being quantified)"
echo "   ${YELLOW}📦 Converted and not quantified:${RESET}   $CONVERTED_NOT_QUANTIFIED (needs quantification)"
echo "   ${YELLOW}🔄 Converting (has SRA):${RESET}           $CONVERTING (SRA → FASTQ in progress)"
echo "   ${YELLOW}⬇️  Downloading (in progress):${RESET}     $DOWNLOADING (downloading SRA)"
echo "   📂 Directories created:            $TOTAL_DIRS"
echo "   ⏸️  Not started:                   $NOT_STARTED"
echo ""

# Verbose mode: show sample lists
if [ "$VERBOSE" -eq 1 ]; then
    echo "   ${GREEN}Completed samples:${RESET}"
    if [ ${#COMPLETED_SAMPLES[@]} -gt 0 ]; then
        printf "     %s\n" "${COMPLETED_SAMPLES[@]}" | head -20
        if [ ${#COMPLETED_SAMPLES[@]} -gt 20 ]; then
            echo "     ... and $(( ${#COMPLETED_SAMPLES[@]} - 20 )) more"
        fi
    else
        echo "     (none)"
    fi
    echo ""
    echo "   ${BLUE}Quantifying samples:${RESET}"
    if [ ${#QUANTIFYING_SAMPLES[@]} -gt 0 ]; then
        printf "     %s\n" "${QUANTIFYING_SAMPLES[@]}"
    else
        echo "     (none)"
    fi
    echo ""
    echo "   ${YELLOW}Not quantified samples:${RESET}"
    if [ ${#NOT_QUANTIFIED_SAMPLES[@]} -gt 0 ]; then
        printf "     %s\n" "${NOT_QUANTIFIED_SAMPLES[@]}" | head -10
        if [ ${#NOT_QUANTIFIED_SAMPLES[@]} -gt 10 ]; then
            echo "     ... and $(( ${#NOT_QUANTIFIED_SAMPLES[@]} - 10 )) more"
        fi
    else
        echo "     (none)"
    fi
    echo ""
fi

# Health warnings section
WARNINGS=0
echo "⚠️  Health Warnings:"
if [ ${#STUCK_SAMPLES[@]} -gt 0 ]; then
    WARNINGS=$((WARNINGS + 1))
    echo "   ${RED}⚠️  Stuck samples detected: ${#STUCK_SAMPLES[@]}${RESET}"
    for stuck_sample in "${!STUCK_SAMPLES[@]}"; do
        stuck_type="${STUCK_SAMPLES[$stuck_sample]}"
        last_activity="${SAMPLE_LAST_ACTIVITY[$stuck_sample]}"
        stuck_minutes=0
        if [ -n "$last_activity" ] && [ "$last_activity" != "0" ]; then
            stuck_minutes=$(( (CURRENT_TIME - last_activity) / 60 ))
        fi
        
        # Get file information for stuck sample
        sample_dir="$FASTQ_DIR/$stuck_sample"
        file_info=""
        if [ -d "$sample_dir" ]; then
            # Check for SRA files
            sra_files=$(find "$sample_dir" -name "*.sra" -type f 2>/dev/null)
            if [ -n "$sra_files" ]; then
                sra_size=$(du -sh $(echo "$sra_files" | head -1) 2>/dev/null | awk '{print $1}')
                sra_count=$(echo "$sra_files" | wc -l)
                file_info=" | SRA: ${sra_count} file(s), ${sra_size}"
                
                # Check if SRA file is growing (recent modification)
                latest_sra=$(echo "$sra_files" | xargs ls -t 2>/dev/null | head -1)
                if [ -n "$latest_sra" ] && [ -f "$latest_sra" ]; then
                    sra_mod_time=$(stat -c %Y "$latest_sra" 2>/dev/null || stat -f %m "$latest_sra" 2>/dev/null || echo "0")
                    if [ "$sra_mod_time" != "0" ]; then
                        sra_age_minutes=$(( (CURRENT_TIME - sra_mod_time) / 60 ))
                        if [ "$sra_age_minutes" -lt 5 ]; then
                            file_info="${file_info} (growing)"
                        fi
                    fi
                fi
            fi
            
            # Check for FASTQ files
            fastq_files=$(find "$sample_dir" -name "*.fastq*" ! -name "*.sra" -type f 2>/dev/null)
            if [ -n "$fastq_files" ]; then
                fastq_count=$(echo "$fastq_files" | wc -l)
                # Calculate total FASTQ size
                fastq_size_total=0
                for fq_file in $fastq_files; do
                    if [ -f "$fq_file" ]; then
                        file_size=$(stat -c %s "$fq_file" 2>/dev/null || stat -f %z "$fq_file" 2>/dev/null || echo "0")
                        fastq_size_total=$((fastq_size_total + file_size))
                    fi
                done
                # Convert to human readable
                if [ "$fastq_size_total" -gt 0 ]; then
                    if [ "$fastq_size_total" -lt 1024 ]; then
                        fastq_size_str="${fastq_size_total}B"
                    elif [ "$fastq_size_total" -lt 1048576 ]; then
                        fastq_size_str="$((fastq_size_total / 1024))K"
                    elif [ "$fastq_size_total" -lt 1073741824 ]; then
                        fastq_size_str="$((fastq_size_total / 1048576))M"
                    else
                        fastq_size_str="$((fastq_size_total / 1073741824))G"
                    fi
                else
                    fastq_size_str="0"
                fi
                if [ -z "$file_info" ]; then
                    file_info=" | FASTQ: ${fastq_count} file(s), ${fastq_size_str}"
                else
                    file_info="${file_info}, FASTQ: ${fastq_count} file(s), ${fastq_size_str}"
                fi
            fi
            
            # Check directory size
            dir_size=$(du -sh "$sample_dir" 2>/dev/null | awk '{print $1}')
            if [ -n "$dir_size" ]; then
                if [ -z "$file_info" ]; then
                    file_info=" | Total: ${dir_size}"
                else
                    file_info="${file_info} | Total: ${dir_size}"
                fi
            fi
        fi
        
        if [ "$stuck_type" = "converting_stuck" ]; then
            echo "     ${RED}🔄 $stuck_sample:${RESET} Stuck converting (SRA present, no FASTQ) for ${stuck_minutes} minutes${file_info}"
        else
            echo "     ${RED}⬇️  $stuck_sample:${RESET} Stuck downloading for ${stuck_minutes} minutes${file_info}"
        fi
    done
    echo ""
fi

# Check disk space
if command -v df >/dev/null 2>&1; then
    DISK_FREE=$(df -BG "$FASTQ_DIR" 2>/dev/null | tail -1 | awk '{print $4}' | sed 's/G//')
    if [ -n "$DISK_FREE" ] && [ "$DISK_FREE" -lt 10 ]; then
        WARNINGS=$((WARNINGS + 1))
        echo "   ${YELLOW}⚠️  Low disk space:${RESET} ${DISK_FREE}GB free (recommended: >10GB)"
        echo ""
    fi
fi

# Check process count
TOTAL_PROCESSES=$((GETFASTQ_COUNT + QUANT_COUNT_PROC))
if [ "$TOTAL_PROCESSES" -gt 50 ]; then
    WARNINGS=$((WARNINGS + 1))
    echo "   ${YELLOW}⚠️  High process count:${RESET} $TOTAL_PROCESSES active processes (may indicate resource exhaustion)"
    echo ""
fi

if [ "$WARNINGS" -eq 0 ]; then
    echo "   ${GREEN}✅ No warnings${RESET}"
    echo ""
fi

# Disk usage breakdown
echo "💾 Disk Usage Breakdown:"
if [ -d "$FASTQ_DIR" ]; then
    # Calculate sizes for different file types
    SRA_SIZE=0
    FASTQ_SIZE=0
    if command -v du >/dev/null 2>&1; then
        # Get size of SRA files
        SRA_FILES=$(find "$FASTQ_DIR" -name "*.sra" -type f 2>/dev/null)
        if [ -n "$SRA_FILES" ]; then
            SRA_SIZE=$(du -sb $(echo "$SRA_FILES" | tr '\n' ' ') 2>/dev/null | awk '{sum+=$1} END {print sum}' || echo "0")
        fi
        
        # Get size of FASTQ files
        FASTQ_FILES=$(find "$FASTQ_DIR" -name "*.fastq*" ! -name "*.sra" -type f 2>/dev/null)
        if [ -n "$FASTQ_FILES" ]; then
            FASTQ_SIZE=$(du -sb $(echo "$FASTQ_FILES" | tr '\n' ' ') 2>/dev/null | awk '{sum+=$1} END {print sum}' || echo "0")
        fi
    fi
    
    TOTAL_FASTQ_SIZE=$((SRA_SIZE + FASTQ_SIZE))
    # Use awk for calculation if bc is not available
    if command -v bc >/dev/null 2>&1; then
        TOTAL_FASTQ_SIZE_GB=$(echo "scale=2; $TOTAL_FASTQ_SIZE / 1073741824" | bc 2>/dev/null || echo "0")
        SRA_SIZE_GB=$(echo "scale=2; $SRA_SIZE / 1073741824" | bc 2>/dev/null || echo "0")
        FASTQ_SIZE_GB=$(echo "scale=2; $FASTQ_SIZE / 1073741824" | bc 2>/dev/null || echo "0")
    else
        TOTAL_FASTQ_SIZE_GB=$(awk "BEGIN {printf \"%.2f\", $TOTAL_FASTQ_SIZE / 1073741824}" 2>/dev/null || echo "0")
        SRA_SIZE_GB=$(awk "BEGIN {printf \"%.2f\", $SRA_SIZE / 1073741824}" 2>/dev/null || echo "0")
        FASTQ_SIZE_GB=$(awk "BEGIN {printf \"%.2f\", $FASTQ_SIZE / 1073741824}" 2>/dev/null || echo "0")
    fi
    
    echo "   FASTQ directory total:         ${TOTAL_FASTQ_SIZE_GB}GB"
    echo "     └─ SRA files:                ${SRA_SIZE_GB}GB"
    echo "     └─ FASTQ files:              ${FASTQ_SIZE_GB}GB"
fi

if [ -d "$QUANT_DIR" ]; then
    QUANT_SIZE=$(du -sb "$QUANT_DIR" 2>/dev/null | awk '{print $1}' || echo "0")
    if command -v bc >/dev/null 2>&1; then
        QUANT_SIZE_GB=$(echo "scale=2; $QUANT_SIZE / 1073741824" | bc 2>/dev/null || echo "0")
    else
        QUANT_SIZE_GB=$(awk "BEGIN {printf \"%.2f\", $QUANT_SIZE / 1073741824}" 2>/dev/null || echo "0")
    fi
    echo "   Quantification results:         ${QUANT_SIZE_GB}GB"
fi
echo ""

# Count files
SRA_COUNT=$(find "$FASTQ_DIR" -name "*.sra" 2>/dev/null | wc -l)
FASTQ_COUNT=$(find "$FASTQ_DIR" -name "*.fastq*" ! -name "*.sra" 2>/dev/null | wc -l)
QUANT_COUNT=$(find "$QUANT_DIR" -name "abundance.tsv" 2>/dev/null | wc -l)
echo "📦 File Counts:"
echo "   SRA files:        $SRA_COUNT"
echo "   FASTQ files:      $FASTQ_COUNT"
echo "   Quantified:      $QUANT_COUNT"
echo ""

# Check running processes with deduplication and sample extraction
echo "🔄 Running Processes:"
GETFASTQ_COUNT=$(ps aux | grep -E "amalgkit getfastq" | grep -v grep | awk '{pid=$2; cmd=""; for(i=11;i<=NF;i++) cmd=cmd" "$i; print pid, cmd}' | sort -u -k1,1 | wc -l)
QUANT_COUNT_PROC=$(ps aux | grep -E "amalgkit quant" | grep -v grep | wc -l)

if [ "$GETFASTQ_COUNT" -gt 0 ] || [ "$QUANT_COUNT_PROC" -gt 0 ]; then
    echo "   Download workers (getfastq):    $GETFASTQ_COUNT"
    if [ "$GETFASTQ_COUNT" -gt 0 ]; then
        # Show processes with sample IDs
        ps aux | grep -E "amalgkit getfastq" | grep -v grep | awk '{pid=$2; cmd=""; for(i=11;i<=NF;i++) cmd=cmd" "$i; print pid, cmd}' | sort -u -k1,1 | head -5 | while read pid cmd; do
            sample=$(echo "$cmd" | grep -o "metadata\.download\.SRR[0-9]*" | sed 's/metadata\.download\.//' | head -1)
            runtime=$(ps -o etime= -p "$pid" 2>/dev/null | tr -d ' ' || echo "N/A")
            if [ -n "$sample" ]; then
                echo "     PID: $pid | Sample: ${BLUE}$sample${RESET} | Runtime: $runtime"
            else
                echo "     PID: $pid | Runtime: $runtime"
            fi
        done
        if [ "$GETFASTQ_COUNT" -gt 5 ]; then
            echo "     ... and $((GETFASTQ_COUNT - 5)) more"
        fi
    fi
    echo "   Quantification workers (quant): $QUANT_COUNT_PROC"
    if [ "$QUANT_COUNT_PROC" -gt 0 ]; then
        ps aux | grep -E "amalgkit quant" | grep -v grep | head -3 | awk '{pid=$2; runtime=$10; sample=""; for(i=11;i<=NF;i++) {if ($i ~ /SRR[0-9]+/) {sample=$i; break}}; if (sample) print "     PID:", pid, "| Sample:", sample, "| Runtime:", runtime; else print "     PID:", pid, "| Runtime:", runtime}'
        if [ "$QUANT_COUNT_PROC" -gt 3 ]; then
            echo "     ... and $((QUANT_COUNT_PROC - 3)) more"
        fi
    fi
else
    echo "   No active amalgkit processes"
fi
echo ""

# Enhanced log activity section with sample IDs and timestamps
if [ -d "$LOG_DIR" ]; then
    echo "📝 Recent Activity (last 2 hours):"
    
    # Getfastq logs - find most recent that was modified in last 2 hours
    MOST_RECENT_GETFASTQ=""
    for log_file in $(ls -t "$LOG_DIR"/*getfastq*.stdout.log "$LOG_DIR"/*getfastq*.log 2>/dev/null | head -5); do
        if [ -f "$log_file" ]; then
            MOD_TIME=$(stat -c %Y "$log_file" 2>/dev/null || stat -f %m "$log_file" 2>/dev/null || echo "0")
            if [ "$MOD_TIME" != "0" ]; then
                TIME_AGO=$((CURRENT_TIME - MOD_TIME))
                # Only use logs from last 2 hours (7200 seconds)
                if [ "$TIME_AGO" -lt 7200 ]; then
                    MOST_RECENT_GETFASTQ="$log_file"
                    break
                fi
            fi
        fi
    done
    
    if [ -n "$MOST_RECENT_GETFASTQ" ] && [ -f "$MOST_RECENT_GETFASTQ" ]; then
        sample=$(basename "$MOST_RECENT_GETFASTQ" | grep -o "SRR[0-9]*" | head -1)
        MOD_TIME=$(stat -c %Y "$MOST_RECENT_GETFASTQ" 2>/dev/null || stat -f %m "$MOST_RECENT_GETFASTQ" 2>/dev/null || echo "0")
        if [ "$MOD_TIME" != "0" ]; then
            TIME_AGO=$((CURRENT_TIME - MOD_TIME))
            if [ "$TIME_AGO" -lt 60 ]; then
                TIME_STR="${TIME_AGO}s ago"
            elif [ "$TIME_AGO" -lt 3600 ]; then
                TIME_STR="$((TIME_AGO / 60))m ago"
            else
                TIME_STR="$((TIME_AGO / 3600))h ago"
            fi
        else
            TIME_STR="unknown"
        fi
        
        echo "   Download (getfastq):"
        if [ -n "$sample" ]; then
            echo "     Sample: ${BLUE}$sample${RESET} | Last update: $TIME_STR"
        fi
        
        # Show more lines, prioritizing errors and important messages
        # First, check for errors
        error_lines=$(tail -50 "$MOST_RECENT_GETFASTQ" 2>/dev/null | grep -iE "(error|failed|exception|traceback|fatal)" | tail -5)
        if [ -n "$error_lines" ]; then
            echo "     ${RED}Errors:${RESET}"
            echo "$error_lines" | sed 's/^/       /'
        fi
        
        # Then show important status messages
        important_lines=$(tail -50 "$MOST_RECENT_GETFASTQ" 2>/dev/null | grep -E "(completed|downloading|converting|end|start|Time elapsed|progress|percent)" | tail -5)
        if [ -n "$important_lines" ]; then
            echo "     Status:"
            echo "$important_lines" | sed 's/^/       /'
        fi
        
        # If no filtered lines, show last 8 lines
        if [ -z "$error_lines" ] && [ -z "$important_lines" ]; then
            tail -8 "$MOST_RECENT_GETFASTQ" 2>/dev/null | sed 's/^/     /' || echo "     (log file empty or unreadable)"
        fi
    else
        echo "   Download (getfastq): (no log files found)"
    fi
    
    # Quant logs - find most recent that was modified in last 2 hours
    MOST_RECENT_QUANT=""
    for log_file in $(ls -t "$LOG_DIR"/*quant*.stdout.log "$LOG_DIR"/*quant*.log 2>/dev/null | head -5); do
        if [ -f "$log_file" ]; then
            MOD_TIME=$(stat -c %Y "$log_file" 2>/dev/null || stat -f %m "$log_file" 2>/dev/null || echo "0")
            if [ "$MOD_TIME" != "0" ]; then
                TIME_AGO=$((CURRENT_TIME - MOD_TIME))
                # Only use logs from last 2 hours (7200 seconds)
                if [ "$TIME_AGO" -lt 7200 ]; then
                    MOST_RECENT_QUANT="$log_file"
                    break
                fi
            fi
        fi
    done
    
    if [ -n "$MOST_RECENT_QUANT" ] && [ -f "$MOST_RECENT_QUANT" ]; then
        sample=$(basename "$MOST_RECENT_QUANT" | grep -o "SRR[0-9]*" | head -1)
        MOD_TIME=$(stat -c %Y "$MOST_RECENT_QUANT" 2>/dev/null || stat -f %m "$MOST_RECENT_QUANT" 2>/dev/null || echo "0")
        if [ "$MOD_TIME" != "0" ]; then
            TIME_AGO=$((CURRENT_TIME - MOD_TIME))
            if [ "$TIME_AGO" -lt 60 ]; then
                TIME_STR="${TIME_AGO}s ago"
            elif [ "$TIME_AGO" -lt 3600 ]; then
                TIME_STR="$((TIME_AGO / 60))m ago"
            else
                TIME_STR="$((TIME_AGO / 3600))h ago"
            fi
        else
            TIME_STR="unknown"
        fi
        
        echo "   Quantification (quant):"
        if [ -n "$sample" ]; then
            echo "     Sample: ${BLUE}$sample${RESET} | Last update: $TIME_STR"
        fi
        
        # Show more lines, prioritizing errors and important messages
        # First, check for errors
        error_lines=$(tail -50 "$MOST_RECENT_QUANT" 2>/dev/null | grep -iE "(error|failed|exception|traceback|fatal|not found|exiting)" | tail -5)
        if [ -n "$error_lines" ]; then
            echo "     ${RED}Errors:${RESET}"
            echo "$error_lines" | sed 's/^/       /'
        fi
        
        # Then show important status messages
        important_lines=$(tail -50 "$MOST_RECENT_QUANT" 2>/dev/null | grep -E "(completed|quantifying|end|start|abundance|output|file)" | tail -5)
        if [ -n "$important_lines" ]; then
            echo "     Status:"
            echo "$important_lines" | sed 's/^/       /'
        fi
        
        # If no filtered lines, show last 8 lines
        if [ -z "$error_lines" ] && [ -z "$important_lines" ]; then
            tail -8 "$MOST_RECENT_QUANT" 2>/dev/null | sed 's/^/     /' || echo "     (log file empty or unreadable)"
        fi
    else
        echo "   Quantification (quant): (no log files found)"
    fi
fi
echo ""

# Progress metrics and time estimates
if [ -d "$LOG_DIR" ] && [ "$TOTAL_SAMPLES" -gt 0 ]; then
    # Try to calculate processing rate from log timestamps
    # Find oldest and newest log files
    OLDEST_LOG=$(find "$LOG_DIR" -name "*.log" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | head -1 | cut -d' ' -f2-)
    NEWEST_LOG=$(find "$LOG_DIR" -name "*.log" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2-)
    
    if [ -n "$OLDEST_LOG" ] && [ -n "$NEWEST_LOG" ] && [ "$OLDEST_LOG" != "$NEWEST_LOG" ]; then
        OLDEST_TIME=$(stat -c %Y "$OLDEST_LOG" 2>/dev/null || stat -f %m "$OLDEST_LOG" 2>/dev/null || echo "0")
        NEWEST_TIME=$(stat -c %Y "$NEWEST_LOG" 2>/dev/null || stat -f %m "$NEWEST_LOG" 2>/dev/null || echo "0")
        
        if [ "$OLDEST_TIME" != "0" ] && [ "$NEWEST_TIME" != "0" ] && [ "$NEWEST_TIME" -gt "$OLDEST_TIME" ]; then
            ELAPSED_SEC=$((NEWEST_TIME - OLDEST_TIME))
            ELAPSED_HOURS=$(echo "scale=2; $ELAPSED_SEC / 3600" | bc 2>/dev/null || echo "0")
            
            # Calculate using awk if bc not available
            if command -v bc >/dev/null 2>&1; then
                ELAPSED_CHECK=$(echo "$ELAPSED_HOURS > 0" | bc 2>/dev/null || echo "0")
            else
                ELAPSED_CHECK=$(awk "BEGIN {print ($ELAPSED_HOURS > 0) ? 1 : 0}" 2>/dev/null || echo "0")
            fi
            
            if [ "$ELAPSED_CHECK" = "1" ] || [ "$(echo "$ELAPSED_HOURS" | awk '{if ($1 > 0) print 1; else print 0}')" = "1" ]; then
                if command -v bc >/dev/null 2>&1; then
                    SAMPLES_PER_HOUR=$(echo "scale=1; $STARTED / $ELAPSED_HOURS" | bc 2>/dev/null || echo "0")
                else
                    SAMPLES_PER_HOUR=$(awk "BEGIN {printf \"%.1f\", $STARTED / $ELAPSED_HOURS}" 2>/dev/null || echo "0")
                fi
                REMAINING=$((TOTAL_SAMPLES - CONVERTED_QUANTIFIED_DELETED))
                
                if command -v bc >/dev/null 2>&1; then
                    RATE_CHECK=$(echo "$SAMPLES_PER_HOUR > 0" | bc 2>/dev/null || echo "0")
                else
                    RATE_CHECK=$(awk "BEGIN {print ($SAMPLES_PER_HOUR > 0) ? 1 : 0}" 2>/dev/null || echo "0")
                fi
                
                if [ "$RATE_CHECK" = "1" ] && [ "$REMAINING" -gt 0 ]; then
                    if command -v bc >/dev/null 2>&1; then
                        HOURS_REMAINING=$(echo "scale=1; $REMAINING / $SAMPLES_PER_HOUR" | bc 2>/dev/null || echo "0")
                    else
                        HOURS_REMAINING=$(awk "BEGIN {printf \"%.1f\", $REMAINING / $SAMPLES_PER_HOUR}" 2>/dev/null || echo "0")
                    fi
                    echo "⏱️  Processing Metrics:"
                    echo "   Processing rate:          ${SAMPLES_PER_HOUR} samples/hour"
                    echo "   Estimated time remaining: ${HOURS_REMAINING} hours"
                    echo ""
                fi
            fi
        fi
    fi
fi

echo "═══════════════════════════════════════════════════════════════════════════════"
