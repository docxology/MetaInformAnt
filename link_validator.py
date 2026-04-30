#!/usr/bin/env python3
"""
Comprehensive Markdown Link Validator

Scans all markdown files in a repository, extracts all links, and validates:
- Internal links: file existence, case-sensitive paths, anchor references
- External links: HTTP status (404 detection)

Generates a detailed report sorted by severity.
"""

import os
import re
import sys
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict
from urllib.parse import urlparse, unquote
from typing import Dict, List, Tuple, Optional, Set
import json

# Regex patterns for markdown links
# [text](url) style
INLINE_LINK_PATTERN = re.compile(r'\[([^\]]*)\]\(([^)]+)\)')
# Reference-style links: [text][ref]
REFERENCE_LINK_PATTERN = re.compile(r'\[([^\]]*)\]\[([^\]]+)\]')
# Reference definitions: [ref]: url
REFERENCE_DEF_PATTERN = re.compile(r'^\[([^\]]+)\]:\s*(.+)$', re.MULTILINE)
# Anchor-only links: [text](#anchor)
ANCHOR_ONLY_PATTERN = re.compile(r'\[([^\]]*)\]\(#([^)]+)\)')

class LinkValidator:
    def __init__(self, repo_path: str, check_external: bool = True, external_timeout: int = 10):
        self.repo_path = Path(repo_path).resolve()
        self.check_external = check_external
        self.external_timeout = external_timeout
        
        # Results storage
        self.all_files: List[Path] = []
        self.file_links: Dict[str, List[Dict]] = defaultdict(list)
        self.broken_internal: List[Dict] = []
        self.broken_external: List[Dict] = []
        self.stats = {
            'total_files': 0,
            'total_links': 0,
            'internal_links': 0,
            'external_links': 0,
            'broken_internal': 0,
            'broken_external': 0,
            'files_with_broken_links': 0
        }
        
        # Reference definitions map (ref -> url)
        self.reference_defs: Dict[str, str] = {}
    
    def find_markdown_files(self) -> List[Path]:
        """Recursively find all .md files in the repository."""
        print(f"Scanning for markdown files in: {self.repo_path}")
        md_files = list(self.repo_path.rglob("*.md"))
        print(f"Found {len(md_files)} markdown files")
        return md_files
    
    def extract_links_from_file(self, file_path: Path) -> List[Dict]:
        """Extract all links from a markdown file."""
        try:
            content = file_path.read_text(encoding='utf-8', errors='replace')
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
            return []
        
        links = []
        relative_to = file_path.parent
        
        # First, collect reference definitions
        self.reference_defs.clear()
        for match in REFERENCE_DEF_PATTERN.finditer(content):
            ref_name = match.group(1)
            ref_url = match.group(2).strip()
            self.reference_defs[ref_name] = ref_url
        
        # Process inline links: [text](url)
        for match in INLINE_LINK_PATTERN.finditer(content):
            text = match.group(1)
            url = match.group(2).strip()
            line_num = content[:match.start()].count('\n') + 1
            
            link_info = {
                'text': text,
                'url': url,
                'line': line_num,
                'type': 'inline',
                'file': str(file_path.relative_to(self.repo_path))
            }
            links.append(link_info)
        
        # Process reference-style links: [text][ref]
        for match in REFERENCE_LINK_PATTERN.finditer(content):
            text = match.group(1)
            ref_name = match.group(2)
            line_num = content[:match.start()].count('\n') + 1
            
            # Resolve reference
            if ref_name in self.reference_defs:
                url = self.reference_defs[ref_name]
            else:
                url = None  # Undefined reference
            
            link_info = {
                'text': text,
                'url': url,
                'ref': ref_name,
                'line': line_num,
                'type': 'reference',
                'file': str(file_path.relative_to(self.repo_path))
            }
            links.append(link_info)
        
        return links
    
    def classify_link(self, url: Optional[str], source_file: Path) -> Tuple[str, Optional[str], Optional[str]]:
        """
        Classify a link and resolve its target.
        
        Returns:
            (category, resolved_path, anchor)
            category: 'internal', 'external', 'anchor', 'undefined', 'mailto', etc.
            resolved_path: absolute Path for internal files, None otherwise
            anchor: anchor string without '#', or None
        """
        if url is None:
            return ('undefined', None, None)
        
        url = url.strip()
        
        # Check for mailto, tel, etc.
        if url.startswith(('mailto:', 'tel:', 'ftp:', 'file:')):
            return ('special', None, None)
        
        # Check for absolute URLs (http, https)
        if url.startswith(('http://', 'https://')):
            return ('external', None, None)
        
        # Anchor-only link: #anchor
        if url.startswith('#'):
            anchor = url[1:]
            return ('anchor', None, anchor)
        
        # Internal link: could be relative path or absolute /docs/ style
        # Remove anchor if present
        parts = url.split('#', 1)
        path_part = parts[0]
        anchor = parts[1] if len(parts) > 1 else None
        
        # Resolve the path
        if url.startswith('/'):
            # Absolute from repo root
            resolved = self.repo_path / path_part.lstrip('/')
        else:
            # Relative to source file
            resolved = (source_file.parent / path_part).resolve()
        
        # Check if it's a directory (likely an index.md)
        if resolved.is_dir():
            # Try index.md
            if (resolved / 'index.md').exists():
                resolved = resolved / 'index.md'
            # Or README.md
            elif (resolved / 'README.md').exists():
                resolved = resolved / 'README.md'
        
        return ('internal', str(resolved), anchor)
    
    def normalize_path(self, path: str) -> str:
        """Normalize path for comparison (resolve symlinks, absolute)."""
        try:
            p = Path(path).resolve()
            return str(p)
        except Exception:
            return str(Path(path).absolute())
    
    def file_exists(self, path: str) -> bool:
        """Check if a file exists."""
        try:
            return Path(path).is_file()
        except Exception:
            return False
    
    def extract_anchors_from_file(self, file_path: str) -> Set[str]:
        """Extract all heading anchors (IDs) from a markdown file."""
        try:
            content = Path(file_path).read_text(encoding='utf-8', errors='replace')
        except Exception:
            return set()
        
        anchors = set()
        
        # Standard markdown headings: ## Heading {#custom-id} or automatic IDs
        # ATX headings: # Heading
        for match in re.finditer(r'^#{1,6}\s+(.+)$', content, re.MULTILINE):
            heading_text = match.group(1).strip()
            # Check for explicit ID: {#my-id}
            id_match = re.search(r'\{#([^}]+)\}', heading_text)
            if id_match:
                anchors.add(id_match.group(1))
            else:
                # Generate automatic ID (lowercase, spaces to hyphens, strip punctuation)
                auto_id = heading_text.lower()
                auto_id = re.sub(r'[^\w\s-]', '', auto_id)
                auto_id = re.sub(r'\s+', '-', auto_id)
                auto_id = auto_id.strip('-')
                if auto_id:
                    anchors.add(auto_id)
        
        # Also check for explicit HTML anchor tags
        for match in re.finditer(r'<a\s+[^>]*id=["\']([^"\']+)["\']', content, re.IGNORECASE):
            anchors.add(match.group(1))
        
        return anchors
    
    def check_external_link(self, url: str) -> bool:
        """Check if an external URL returns a non-404 status."""
        # Use curl to check HEAD first, fallback to GET
        try:
            # First try HEAD
            result = subprocess.run(
                ['curl', '-s', '-o', '/dev/null', '-w', '%{http_code}', 
                 '--max-time', str(self.external_timeout), '-I', url],
                capture_output=True, text=True, timeout=self.external_timeout + 5
            )
            status = result.stdout.strip()
            if status and status not in ('000', '400', '401', '403', '404', '500', '502', '503'):
                return True
            
            # If HEAD fails or shows error, try GET with range
            result = subprocess.run(
                ['curl', '-s', '-o', '/dev/null', '-w', '%{http_code}',
                 '--max-time', str(self.external_timeout), '--range', '0-0', url],
                capture_output=True, text=True, timeout=self.external_timeout + 5
            )
            status = result.stdout.strip()
            return status and status not in ('400', '401', '403', '404', '500', '502', '503')
        except subprocess.TimeoutExpired:
            return False
        except Exception as e:
            print(f"Error checking external link {url}: {e}")
            return False
    
    def validate_all(self):
        """Main validation loop."""
        self.all_files = self.find_markdown_files()
        self.stats['total_files'] = len(self.all_files)
        
        print("\nExtracting links from all files...")
        for i, file_path in enumerate(self.all_files, 1):
            if i % 100 == 0:
                print(f"  Processed {i}/{len(self.all_files)} files...")
            
            links = self.extract_links_from_file(file_path)
            self.file_links[str(file_path.relative_to(self.repo_path))] = links
            self.stats['total_links'] += len(links)
        
        print(f"\nTotal links extracted: {self.stats['total_links']}")
        print("Validating links...")
        
        # Cache for file anchors to avoid repeated extraction
        anchors_cache = {}
        
        # Validate each file's links
        processed = 0
        for file_path_str, links in self.file_links.items():
            processed += 1
            if processed % 100 == 0:
                print(f"  Validated {processed}/{len(self.file_links)} files...")
            
            source_file = self.repo_path / file_path_str
            
            for link in links:
                url = link.get('url')
                if url is None:
                    # Undefined reference - broken
                    self.broken_internal.append({
                        'file': file_path_str,
                        'line': link['line'],
                        'text': link.get('text', ''),
                        'url': link.get('ref', '(undefined reference)'),
                        'severity': 'HIGH',
                        'reason': 'Undefined reference'
                    })
                    self.stats['broken_internal'] += 1
                    continue
                
                # Categorize link
                category, resolved_path, anchor = self.classify_link(url, source_file)
                
                if category == 'internal':
                    self.stats['internal_links'] += 1
                    # Check if target file exists
                    if not self.file_exists(resolved_path):
                        self.broken_internal.append({
                            'file': file_path_str,
                            'line': link['line'],
                            'text': link.get('text', ''),
                            'url': url,
                            'target': resolved_path,
                            'severity': 'HIGH',
                            'reason': 'Target file does not exist'
                        })
                        self.stats['broken_internal'] += 1
                    elif anchor:
                        # Check anchor exists in target file
                        target_abs = resolved_path
                        if target_abs not in anchors_cache:
                            anchors_cache[target_abs] = self.extract_anchors_from_file(target_abs)
                        
                        if anchor not in anchors_cache[target_abs]:
                            self.broken_internal.append({
                                'file': file_path_str,
                                'line': link['line'],
                                'text': link.get('text', ''),
                                'url': url,
                                'target': resolved_path,
                                'anchor': anchor,
                                'severity': 'MEDIUM',
                                'reason': f'Anchor "{anchor}" not found in target file'
                            })
                            self.stats['broken_internal'] += 1
                
                elif category == 'external' and self.check_external:
                    self.stats['external_links'] += 1
                    if not self.check_external_link(url):
                        self.broken_external.append({
                            'file': file_path_str,
                            'line': link['line'],
                            'text': link.get('text', ''),
                            'url': url,
                            'severity': 'MEDIUM',
                            'reason': 'External URL returns 404 or unreachable'
                        })
                        self.stats['broken_external'] += 1
                
                elif category == 'anchor':
                    # Standalone anchor link - need to check if anchor exists in current file
                    # Get anchors for current file
                    current_file_abs = str(source_file.resolve())
                    if current_file_abs not in anchors_cache:
                        anchors_cache[current_file_abs] = self.extract_anchors_from_file(source_file)
                    
                    if anchor not in anchors_cache[current_file_abs]:
                        self.broken_internal.append({
                            'file': file_path_str,
                            'line': link['line'],
                            'text': link.get('text', ''),
                            'url': url,
                            'anchor': anchor,
                            'severity': 'LOW',
                            'reason': f'Anchor "{anchor}" not found in current file'
                        })
                        self.stats['broken_internal'] += 1
                elif category == 'undefined':
                    self.stats['broken_internal'] += 1
                # Ignore special categories (mailto, etc.)
        
        # Calculate files with broken links
        files_with_broken = set()
        for item in self.broken_internal + self.broken_external:
            files_with_broken.add(item['file'])
        self.stats['files_with_broken_links'] = len(files_with_broken)
    
    def generate_report(self, output_path: Optional[Path] = None):
        """Generate comprehensive broken-link report."""
        report_lines = []
        
        report_lines.append("# Markdown Link Validation Report\n")
        report_lines.append(f"**Repository:** {self.repo_path}\n")
        report_lines.append(f"**Scanned:** {self.stats['total_files']} markdown files\n")
        report_lines.append(f"**Total links found:** {self.stats['total_links']}\n")
        report_lines.append(f"**Internal links:** {self.stats['internal_links']}\n")
        if self.check_external:
            report_lines.append(f"**External links:** {self.stats['external_links']}\n")
        report_lines.append("\n## Summary\n")
        report_lines.append(f"- **Broken internal links:** {self.stats['broken_internal']}")
        report_lines.append(f"- **Broken external links:** {self.stats['broken_external']}")
        report_lines.append(f"- **Files with broken links:** {self.stats['files_with_broken_links']}\n")
        
        # Severity breakdown
        severity_counts = defaultdict(int)
        for item in self.broken_internal + self.broken_external:
            severity_counts[item['severity']] += 1
        
        report_lines.append("### Severity Breakdown\n")
        for severity in ['HIGH', 'MEDIUM', 'LOW']:
            report_lines.append(f"- {severity}: {severity_counts.get(severity, 0)}\n")
        
        # Detailed broken links - sorted by severity then file
        report_lines.append("\n## Broken Links (Sorted by Severity)\n")
        report_lines.append("### HIGH - Missing files or undefined references\n")
        high_items = [item for item in self.broken_internal if item['severity'] == 'HIGH']
        high_items.sort(key=lambda x: (x['file'], x['line']))
        for item in high_items:
            report_lines.append(
                f"- **{item['file']}** (line {item['line']}): "
                f"`{item['url']}` → {item['reason']} "
                f"[text: '{item.get('text', '')}']"
            )
        
        report_lines.append("\n### MEDIUM - Broken anchors or inaccessible external URLs\n")
        med_items = [item for item in self.broken_internal + self.broken_external 
                     if item['severity'] == 'MEDIUM']
        med_items.sort(key=lambda x: (x['file'], x['line']))
        for item in med_items:
            target_info = f" (anchor: {item.get('anchor')})" if item.get('anchor') else ""
            report_lines.append(
                f"- **{item['file']}** (line {item['line']}): "
                f"`{item['url']}` → {item['reason']}{target_info} "
                f"[text: '{item.get('text', '')}']"
            )
        
        report_lines.append("\n### LOW - Minor anchor issues in same file\n")
        low_items = [item for item in self.broken_internal if item['severity'] == 'LOW']
        low_items.sort(key=lambda x: (x['file'], x['line']))
        for item in low_items:
            report_lines.append(
                f"- **{item['file']}** (line {item['line']}): "
                f"`{item['url']}` → {item['reason']} "
                f"[text: '{item.get('text', '')}']"
            )
        
        # Per-file breakdown
        report_lines.append("\n## Per-File Broken Link Count\n")
        per_file = defaultdict(int)
        for item in self.broken_internal + self.broken_external:
            per_file[item['file']] += 1
        
        for file_path, count in sorted(per_file.items(), key=lambda x: (-x[1], x[0])):
            report_lines.append(f"- {file_path}: {count} broken link(s)")
        
        report_text = '\n'.join(report_lines)
        
        if output_path:
            output_path.write_text(report_text, encoding='utf-8')
            print(f"\nReport saved to: {output_path}")
        else:
            print("\n" + report_text)
        
        return report_text
    
    def save_json_report(self, output_path: Path):
        """Save detailed JSON report for programmatic consumption."""
        data = {
            'repository': str(self.repo_path),
            'statistics': self.stats,
            'broken_internal_links': self.broken_internal,
            'broken_external_links': self.broken_external,
            'timestamp': subprocess.getoutput('date -Iseconds')
        }
        output_path.write_text(json.dumps(data, indent=2, ensure_ascii=False))
        print(f"JSON report saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Validate markdown documentation links')
    parser.add_argument('--repo', default='/home/trim/Documents/Git/MetaInformAnt',
                       help='Repository path (default: provided workspace)')
    parser.add_argument('--no-external', action='store_true',
                       help='Skip checking external URLs (faster)')
    parser.add_argument('--timeout', type=int, default=10,
                       help='Timeout for external URL checks in seconds (default: 10)')
    parser.add_argument('--output', default='LINK_VALIDATION_REPORT.md',
                       help='Output markdown report filename')
    parser.add_argument('--json', default='link_validation_report.json',
                       help='Output JSON report filename')
    parser.add_argument('--json-only', action='store_true',
                       help='Only generate JSON report (no console output)')
    
    args = parser.parse_args()
    
    validator = LinkValidator(
        repo_path=args.repo,
        check_external=not args.no_external,
        external_timeout=args.timeout
    )
    
    try:
        validator.validate_all()
        
        if not args.json_only:
            output_path = Path(args.output)
            validator.generate_report(output_path)
        
        json_path = Path(args.json)
        validator.save_json_report(json_path)
        
        print("\n" + "="*60)
        print("VALIDATION COMPLETE")
        print("="*60)
        print(f"Files scanned: {validator.stats['total_files']}")
        print(f"Total links: {validator.stats['total_links']}")
        print(f"Broken internal: {validator.stats['broken_internal']}")
        print(f"Broken external: {validator.stats['broken_external']}")
        print(f"Files affected: {validator.stats['files_with_broken_links']}")
        print()
        
        if validator.stats['broken_internal'] > 0 or validator.stats['broken_external'] > 0:
            print("ISSUES FOUND - See detailed report above or in output files")
            return 1
        else:
            print("SUCCESS - All links validated successfully")
            return 0
            
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        return 130
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
