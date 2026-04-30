#!/usr/bin/env python3
"""
Core Documentation Cross-Check

Compares documented API in docs/core/*.md against implementation in src/metainformant/core/*.py.
Generates a report of mismatches, missing, and extra symbols.
"""

import ast
import re
import sys
from pathlib import Path
from collections import defaultdict

# --- Configuration ---
DOCS_DIR = Path("docs/core")
SRC_DIR = Path("src/metainformant/core")
# Exclude these modules from code inventory (internal/private)
EXCLUDE_CODE_MODULES = {'__init__'}  # skip __init__ files for module assignment
# Skip symbols starting with underscore (private) on both sides
SKIP_PRIVATE = True

# --- Helper functions ---

def is_identifier(s):
    return re.match(r'^[A-Za-z_]\w*$', s) is not None

def clean_whitespace(s):
    return re.sub(r'\s+', ' ', s).strip()

def parse_doc_signature_to_params(sig_str):
    """
    Parse a documented function/constructor signature into a normalized
    parameter structure: {'pos': [(name, default_str)], 'vararg': name|None,
                         'kwonly': [(name, default_str)], 'kwarg': name|None}
    Returns None if parsing fails.
    """
    # Remove return annotation if present
    before = sig_str.split('->')[0].strip()
    # Find content inside outermost parentheses
    start = before.find('(')
    if start == -1:
        # No parentheses -> no params
        return {'pos': [], 'vararg': None, 'kwonly': [], 'kwarg': None}
    depth = 0
    end = -1
    for i in range(start, len(before)):
        c = before[i]
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
            if depth == 0:
                end = i
                break
    if end == -1:
        param_str = before[start+1:]
    else:
        param_str = before[start+1:end].strip()
    # Collapse whitespace
    param_str = re.sub(r'\s+', ' ', param_str)
    if not param_str:
        return {'pos': [], 'vararg': None, 'kwonly': [], 'kwarg': None}
    # Build dummy function
    dummy_src = f"def _({param_str}): pass"
    try:
        tree = ast.parse(dummy_src)
        func = tree.body[0]
        args = func.args
        # Positional args
        pos_args = args.args
        defaults = args.defaults
        num_pos = len(pos_args)
        num_defaults = len(defaults)
        pos_with_defaults = []
        for i, arg in enumerate(pos_args):
            default_idx = i - (num_pos - num_defaults)
            if default_idx >= 0:
                default_node = defaults[default_idx]
                default_str = ast.unparse(default_node)
            else:
                default_str = ''
            pos_with_defaults.append((arg.arg, default_str))
        vararg = args.vararg.arg if args.vararg else None
        kwonly = []
        for i, kwarg in enumerate(args.kwonlyargs):
            default_node = args.kw_defaults[i] if i < len(args.kw_defaults) else None
            default_str = ast.unparse(default_node) if default_node is not None else ''
            kwonly.append((kwarg.arg, default_str))
        kwarg = args.kwarg.arg if args.kwarg else None
        return {'pos': pos_with_defaults, 'vararg': vararg, 'kwonly': kwonly, 'kwarg': kwarg}
    except Exception as e:
        return None

def params_to_tuple(p):
    """Convert params structure to a hashable tuple for comparison."""
    if p is None:
        return None
    pos = tuple((name, default) for name, default in p['pos'])
    vararg = p['vararg'] if p['vararg'] is not None else ''
    kwonly = tuple((name, default) for name, default in p['kwonly'])
    kwarg = p['kwarg'] if p['kwarg'] is not None else ''
    return (pos, vararg, kwonly, kwarg)

def params_equal(p1, p2):
    if p1 is None or p2 is None:
        return p1 == p2  # both None or one None
    return params_to_tuple(p1) == params_to_tuple(p2)

# --- Doc parsing ---

def extract_symbols_from_doc(doc_path):
    symbols = []
    with open(doc_path, encoding='utf-8') as f:
        lines = f.readlines()

    current_class = None
    current_class_level = None
    in_table = False
    table_name_idx = None
    table_sig_idx = None
    last_pipe_line = None   # potential table header line
    in_code_block = False
    code_lines = []
    code_start_line = None

    def process_code_block():
        nonlocal code_lines, code_start_line, current_class
        block_text = ''.join(code_lines)
        # Skip if block contains a function or class definition (likely an example)
        if re.search(r'^\s*(async\s+)?def\s+', block_text, re.MULTILINE) or re.search(r'^\s*class\s+', block_text, re.MULTILINE):
            return
        stripped_lines = [l.rstrip() for l in code_lines if l.strip() != '']
        if not stripped_lines:
            return
        first_line = stripped_lines[0].strip()
        # Check if first line starts with an identifier followed by '(' (constructor-like)
        if re.match(r'^[A-Za-z_]\w*\s*\(', first_line):
            # Possibly a constructor signature
            # Full signature (joined)
            sig = ' '.join(l.strip() for l in stripped_lines)
            name = first_line.split('(')[0].strip()
            # Validate identifier
            if is_identifier(name):
                qualified_name = name  # usually class name
                if not name.startswith('_'):
                    symbols.append({
                        'name': name,
                        'qualified_name': qualified_name,
                        'signature': sig,
                        'params': None,  # will be parsed later
                        'doc_file': str(doc_path),
                        'line': code_start_line,
                    })

    for lineno, line in enumerate(lines, start=1):
        stripped = line.strip()

        # Code fence handling
        if stripped.startswith('```'):
            if not in_code_block:
                in_code_block = True
                code_start_line = lineno
                code_lines = []
            else:
                in_code_block = False
                process_code_block()
            continue

        if in_code_block:
            code_lines.append(line)
            continue

        # Track potential table header lines (non-fence, non-code, starting with | and ending with |, not a separator)
        if (not in_table and stripped.startswith('|') and stripped.endswith('|') and
                '|' in stripped[1:-1] and not re.match(r'^\|[\s:-]+\|$', stripped)):
            last_pipe_line = line
            # Don't continue; let further checks happen? Actually this line might be the header itself. We set last_pipe_line and then if next line is separator we'll handle. We'll continue to next iteration so that we don't treat this as data row yet.
            continue

        # Check for table separator line
        if re.match(r'^\|[\s:-]+\|$', stripped):
            if last_pipe_line is not None:
                # Parse header
                header_line = last_pipe_line
                header_cells = [c.strip() for c in header_line.strip().split('|')[1:-1]]
                name_idx = None
                sig_idx = None
                for idx, cell in enumerate(header_cells):
                    cl = cell.lower()
                    if 'method' in cl or 'name' in cl or 'function' in cl:
                        name_idx = idx
                    if 'signature' in cl:
                        sig_idx = idx
                if name_idx is not None and sig_idx is not None:
                    in_table = True
                    table_name_idx = name_idx
                    table_sig_idx = sig_idx
                # Reset last_pipe_line
                last_pipe_line = None
            continue

        if in_table:
            # Check for row
            if stripped.startswith('|') and stripped.endswith('|') and '|' in stripped[1:-1]:
                cells = [c.strip() for c in stripped.split('|')[1:-1]]
                if len(cells) > max(table_name_idx, table_sig_idx):
                    name_cell = cells[table_name_idx]
                    sig_cell = cells[table_sig_idx]
                    name_match = re.search(r'`([^`]+)`', name_cell)
                    sig_match = re.search(r'`([^`]+)`', sig_cell)
                    if name_match:
                        name = name_match.group(1).strip()
                        # Validate identifier
                        if is_identifier(name):
                            signature = sig_match.group(1).strip() if sig_match else ''
                            # Determine qualified name
                            if current_class:
                                qualified_name = f"{current_class}.{name}"
                            else:
                                qualified_name = name
                            if not name.startswith('_'):
                                symbols.append({
                                    'name': name,
                                    'qualified_name': qualified_name,
                                    'signature': signature,
                                    'params': None,
                                    'doc_file': str(doc_path),
                                    'line': lineno,
                                })
                continue
            else:
                in_table = False
                # Continue to handle this line normally (could be heading)
                # Do not continue, fall through

        # Heading detection
        if line.startswith('#'):
            m = re.match(r'^(#+)\s*(.*)', line)
            if m:
                level = len(m.group(1))
                text = m.group(2).strip()
                bt_match = re.search(r'`([^`]+)`', text)
                if bt_match:
                    symbol_text = bt_match.group(1).strip()
                    # Determine if contains '(' -> function/method/constructor
                    if '(' in symbol_text:
                        base_name = symbol_text.split('(')[0].strip()
                        if is_identifier(base_name):
                            # method or function
                            if current_class is not None and level > current_class_level:
                                qualified_name = f"{current_class}.{base_name}"
                            else:
                                qualified_name = base_name
                            if not base_name.startswith('_'):
                                symbols.append({
                                    'name': base_name,
                                    'qualified_name': qualified_name,
                                    'signature': symbol_text,
                                    'params': None,
                                    'doc_file': str(doc_path),
                                    'line': lineno,
                                })
                    else:
                        # Plain identifier: treat as class name
                        class_name = symbol_text
                        if is_identifier(class_name) and not class_name.startswith('_'):
                            symbols.append({
                                'name': class_name,
                                'qualified_name': class_name,
                                'signature': '',
                                'params': None,
                                'doc_file': str(doc_path),
                                'line': lineno,
                            })
                        # Update class context
                        current_class = class_name
                        current_class_level = level
                else:
                    # No backticks in heading: might be section header
                    # If heading level <= current_class_level, reset class context
                    if current_class is not None and level <= current_class_level:
                        current_class = None
            continue  # processed heading

    return symbols

# --- Code parsing ---

class DefinitionVisitor(ast.NodeVisitor):
    def __init__(self, filename, lines):
        self.filename = filename
        self.lines = lines
        self.symbols = []
        self.class_stack = []
        self.class_inits = {}  # class_name -> (node, lineno)

    def visit_ClassDef(self, node):
        class_name = node.name
        if SKIP_PRIVATE and class_name.startswith('_'):
            return
        # Add class symbol
        self.symbols.append({
            'name': class_name,
            'qualified_name': class_name,
            'type': 'class',
            'signature_raw': '',
            'params': None,
            'file': self.filename,
            'line': node.lineno,
        })
        self.class_stack.append(class_name)
        self.generic_visit(node)
        self.class_stack.pop()

    def visit_FunctionDef(self, node):
        self._handle_function(node)

    def visit_AsyncFunctionDef(self, node):
        self._handle_function(node)

    def _handle_function(self, node):
        name = node.name
        if SKIP_PRIVATE and name.startswith('_'):
            if name == '__init__' and self.class_stack:
                # Store constructor for later
                class_name = self.class_stack[-1]
                self.class_inits[class_name] = (node, self.filename, node.lineno)
            return

        is_method = bool(self.class_stack)
        if is_method:
            qualified_name = f"{self.class_stack[-1]}.{name}"
        else:
            qualified_name = name

        raw_sig = self._get_signature_raw(node)
        params = self._parse_params(node, is_method=is_method)

        self.symbols.append({
            'name': name,
            'qualified_name': qualified_name,
            'type': 'method' if is_method else 'function',
            'signature_raw': raw_sig,
            'params': params,
            'file': self.filename,
            'line': node.lineno,
        })

    def _get_signature_raw(self, node):
        start = node.lineno - 1
        # Check first statement line
        if node.body:
            first_stmt_line = node.body[0].lineno - 1
        else:
            first_stmt_line = start
        if first_stmt_line > start:
            # multi-line signature
            sig_lines = self.lines[start:first_stmt_line]
            return '\n'.join(l.rstrip() for l in sig_lines).strip()
        else:
            # single line; cut at colon
            line = self.lines[start]
            colon_idx = line.find(':')
            if colon_idx != -1:
                return line[:colon_idx+1].strip()
            else:
                return line.strip()

    def _parse_params(self, node, is_method):
        args = node.args
        # Positional args
        pos_args = args.args
        if is_method and pos_args:
            # skip first argument (self/cls)
            pos_args = pos_args[1:]
        num_pos = len(pos_args)
        defaults = args.defaults
        num_defaults = len(defaults)
        pos_with_defaults = []
        for i, arg in enumerate(pos_args):
            default_idx = i - (num_pos - num_defaults)
            if default_idx >= 0:
                default_node = defaults[default_idx]
                try:
                    default_str = ast.unparse(default_node)
                except AttributeError:
                    default_str = ''
            else:
                default_str = ''
            pos_with_defaults.append((arg.arg, default_str))
        # vararg
        vararg = args.vararg.arg if args.vararg else None
        # keyword-only
        kwonly = []
        kwonly_defaults = args.kw_defaults
        for i, kwarg in enumerate(args.kwonlyargs):
            default_node = kwonly_defaults[i] if i < len(kwonly_defaults) else None
            default_str = ast.unparse(default_node) if default_node is not None else ''
            kwonly.append((kwarg.arg, default_str))
        # kwarg
        kwarg = args.kwarg.arg if args.kwarg else None
        return {'pos': pos_with_defaults, 'vararg': vararg, 'kwonly': kwonly, 'kwarg': kwarg}

def extract_code_symbols(src_root):
    all_symbols = []
    py_files = list(src_root.rglob("*.py"))
    for py_file in py_files:
        # Skip files in __pycache__
        if '__pycache__' in py_file.parts:
            continue
        with open(py_file, encoding='utf-8') as f:
            content = f.read()
            lines = content.splitlines()
        try:
            tree = ast.parse(content, filename=str(py_file))
        except SyntaxError as e:
            print(f"Warning: skipping {py_file} due to syntax error: {e}", file=sys.stderr)
            continue
        visitor = DefinitionVisitor(str(py_file), lines)
        visitor.visit(tree)
        # After visit, add constructor symbols
        for class_name, (init_node, fname, lineno) in visitor.class_inits.items():
            # Create constructor symbol
            raw_sig = visitor._get_signature_raw(init_node)
            params = visitor._parse_params(init_node, is_method=True)  # skips self
            all_symbols.append({
                'name': class_name,
                'qualified_name': class_name,
                'type': 'constructor',
                'signature_raw': raw_sig,
                'params': params,
                'file': fname,
                'line': lineno,
            })
        all_symbols.extend(visitor.symbols)
    return all_symbols

# --- Main comparison ---

def main():
    print("=== Core Documentation Cross-Check ===\n")

    # 1. Parse docs
    doc_files = list(DOCS_DIR.glob("*.md"))
    doc_symbols = []
    for doc_file in doc_files:
        syms = extract_symbols_from_doc(doc_file)
        doc_symbols.extend(syms)

    # Filter private symbols from docs if needed
    if SKIP_PRIVATE:
        doc_symbols = [s for s in doc_symbols if not s['name'].startswith('_')]

    # Parse params for doc symbols (those with signature) after extraction, because we need to collapse etc.
    for s in doc_symbols:
        if s['signature'] and s['params'] is None:
            s['params'] = parse_doc_signature_to_params(s['signature'])

    # 2. Parse code
    code_symbols = extract_code_symbols(SRC_DIR)
    # Filter private from code as well (already skipped most, but ensure)
    if SKIP_PRIVATE:
        code_symbols = [s for s in code_symbols if not s['name'].startswith('_')]

    # 3. Build lookup
    code_by_qname = {s['qualified_name']: s for s in code_symbols}
    # Also build simple name mapping for fallback? But we rely on qualified_name primarily.

    # 4. Compare
    matched_qnames = set()
    missing_in_code = []   # doc symbols not found in code
    mismatches = []       # signature mismatches

    for ds in doc_symbols:
        qname = ds['qualified_name']
        if qname in code_by_qname:
            cs = code_by_qname[qname]
            # Compare signatures if both have params
            doc_params = ds.get('params')
            code_params = cs.get('params')
            if doc_params is not None and code_params is not None:
                if params_equal(doc_params, code_params):
                    matched_qnames.add(qname)
                else:
                    mismatches.append({
                        'symbol': qname,
                        'doc_file': ds['doc_file'],
                        'doc_line': ds['line'],
                        'doc_signature': ds['signature'],
                        'code_file': cs['file'],
                        'code_line': cs['line'],
                        'code_signature': cs['signature_raw'],
                    })
            else:
                # If either has no params (e.g., class with empty signature), treat as match
                matched_qnames.add(qname)
        else:
            missing_in_code.append(ds)

    # Missing in docs: code symbols not matched by any doc symbol
    missing_in_docs = [cs for cs in code_symbols if cs['qualified_name'] not in matched_qnames]

    # 5. Per-module accuracy (based on code modules)
    # Module id: stem of the python file, excluding __init__
    module_stats = defaultdict(lambda: {'total': 0, 'covered': 0})
    for cs in code_symbols:
        modname = Path(cs['file']).stem
        if modname in EXCLUDE_CODE_MODULES:
            continue
        module_stats[modname]['total'] += 1
        if cs['qualified_name'] in matched_qnames:
            module_stats[modname]['covered'] += 1

    # 6. Print report
    print(f"Total documented symbols: {len(doc_symbols)}")
    print(f"Total implemented symbols: {len(code_symbols)}")
    print(f"Documented but missing in code: {len(missing_in_code)}")
    print(f"Implemented but undocumented: {len(missing_in_docs)}")
    print(f"Signature mismatches: {len(mismatches)}\n")

    print("--- Per-Module Accuracy (Implementation Coverage) ---")
    for mod, stats in sorted(module_stats.items()):
        total = stats['total']
        if total == 0:
            continue
        covered = stats['covered']
        acc = (covered / total) * 100
        print(f"  {mod:30s}  {covered}/{total}  ({acc:.1f}%)")

    if missing_in_code:
        print("\n### Documented symbols missing in code ###")
        for s in missing_in_code:
            print(f"  {s['qualified_name']}  (in {Path(s['doc_file']).name} line {s['line']})")

    if missing_in_docs:
        print("\n### Implemented symbols missing from documentation ###")
        for s in missing_in_docs:
            print(f"  {s['qualified_name']}  ({s['file']}:{s['line']})")

    if mismatches:
        print("\n### Signature mismatches ###")
        for m in mismatches:
            print(f"  Symbol: {m['symbol']}")
            print(f"    Doc:    {m['doc_signature']}  ({Path(m['doc_file']).name} line {m['doc_line']})")
            print(f"    Code:   {m['code_signature']}  ({Path(m['code_file']).name} line {m['code_line']})")
            print()

if __name__ == '__main__':
    main()
