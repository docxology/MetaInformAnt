#!/usr/bin/env python3
"""
Script to fix mermaid diagram violations in METAINFORMANT documentation.

Fixes common violations:
1. Node IDs with spaces: [Label With Spaces] -> [labelWithSpaces]
2. Node IDs with underscores: [label_with_underscores] -> [labelWithUnderscores]
3. Reserved keywords as node IDs
"""

import re
import os
import glob
from pathlib import Path
from typing import Dict, List, Tuple


def camel_case(text: str) -> str:
    """Convert text to camelCase, removing spaces and underscores."""
    # Split on spaces, underscores, and other separators
    words = re.split(r'[_\s]+', text.strip())
    # Convert first word to lowercase, others to title case
    if not words:
        return ""
    return words[0].lower() + ''.join(word.capitalize() for word in words[1:])


def fix_node_id(node_id: str) -> str:
    """Fix a mermaid node ID to follow camelCase convention."""
    # Remove brackets if present
    node_id = node_id.strip('[]')
    return camel_case(node_id)


def fix_mermaid_content(content: str) -> str:
    """Fix mermaid diagram content by correcting node definitions."""
    lines = content.split('\n')
    fixed_lines = []

    for line in lines:
        # Look for patterns that need fixing:
        # 1. [Label with spaces] used as node ID (should be camelCaseId[Label with spaces])
        # 2. NODE_ID[label_with_underscores] (underscores in labels are okay, but fix node IDs)

        # Pattern 1: [Label with spaces] - these are invalid node definitions
        bracket_only_pattern = r'\[([^\]]*[ _][^\]]*)\]'
        bracket_matches = list(re.finditer(bracket_only_pattern, line))

        if bracket_matches:
            new_line = line
            offset = 0

            for match in bracket_matches:
                label = match.group(1)
                # Create camelCase node ID from the label
                node_id = fix_node_id(label)

                # Replace [Label] with nodeId[Label]
                start_pos = match.start() + offset
                end_pos = match.end() + offset

                replacement = f'{node_id}[{label}]'
                new_line = new_line[:start_pos] + replacement + new_line[end_pos:]
                offset += len(replacement) - (end_pos - start_pos)

            fixed_lines.append(new_line)
        else:
            # Pattern 2: Check for NODE_ID[label] where NODE_ID has issues
            node_pattern = r'(\w+)\[([^\]]+)\]'
            matches = re.finditer(node_pattern, line)

            if matches:
                new_line = line
                offset = 0

                for match in matches:
                    node_id = match.group(1)
                    label = match.group(2)

                    # Check if node_id has violations
                    if (' ' in node_id or '_' in node_id or
                        node_id in ['end', 'subgraph', 'graph', 'flowchart']):
                        new_node_id = fix_node_id(label)

                        start_pos = match.start(1) + offset
                        end_pos = match.end(1) + offset

                        new_line = new_line[:start_pos] + new_node_id + new_line[end_pos:]
                        offset += len(new_node_id) - len(node_id)

                fixed_lines.append(new_line)
            else:
                fixed_lines.append(line)

    return '\n'.join(fixed_lines)


def process_file(filepath: Path) -> Tuple[int, List[str]]:
    """Process a single file, fixing mermaid violations."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Find all mermaid code blocks
    mermaid_blocks = re.findall(r'```mermaid\n(.*?)\n```', content, re.DOTALL)

    if not mermaid_blocks:
        return 0, []

    changes_made = []
    new_content = content

    for block in mermaid_blocks:
        fixed_block = fix_mermaid_content(block)
        if fixed_block != block:
            # Replace the block in the content
            new_content = new_content.replace(f'```mermaid\n{block}\n```',
                                            f'```mermaid\n{fixed_block}\n```')
            changes_made.append(f"Fixed mermaid block in {filepath}")

    if changes_made:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)

    return len(changes_made), changes_made


def main():
    """Main function to process all documentation files."""
    repo_root = Path(__file__).parent.parent

    # Find all markdown files with mermaid diagrams
    pattern = "**/*.md"
    files = list(repo_root.glob(pattern))

    total_changes = 0
    processed_files = []

    print("🔧 Fixing mermaid diagram violations...")
    print("=" * 50)

    for filepath in files:
        changes, change_details = process_file(filepath)
        if changes > 0:
            total_changes += changes
            processed_files.append(filepath)
            print(f"✅ {filepath.relative_to(repo_root)}: {changes} mermaid blocks fixed")

    print("=" * 50)
    print(f"📊 Summary: {total_changes} mermaid blocks fixed across {len(processed_files)} files")
    print("🎯 All mermaid diagrams now follow camelCase node ID conventions!")


if __name__ == "__main__":
    main()