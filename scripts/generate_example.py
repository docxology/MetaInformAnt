#!/usr/bin/env python3
"""Generate new METAINFORMANT examples from templates.

This script creates new example files based on domain-specific templates,
ensuring consistent structure and best practices.

Usage:
    python scripts/generate_example.py <domain> <name> [--description TEXT] [--features LIST]

Arguments:
    domain: Target domain (dna, rna, gwas, etc.)
    name: Example name (will create example_<name>.py)

Options:
    --description: Brief description of what the example demonstrates
    --features: Comma-separated list of features demonstrated
    --template: Specific template to use (default: auto-detect)
    --dry-run: Show what would be created without creating files
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List


class ExampleGenerator:
    """Generate METAINFORMANT examples from templates."""

    def __init__(self):
        self.templates_dir = Path("examples/templates")
        self.examples_dir = Path("examples")

        # Domain-specific configurations
        self.domain_configs = {
            "core": {
                "template": "base_template.py",
                "description": "Core utilities demonstration",
                "features": ["Configuration management", "I/O operations", "Logging setup"]
            },
            "dna": {
                "template": "dna_template.py",
                "description": "DNA sequence analysis",
                "features": ["Sequence processing", "Analysis algorithms", "Result visualization"]
            },
            "rna": {
                "template": "rna_template.py",
                "description": "RNA analysis workflows",
                "features": ["Transcriptomic analysis", "Expression quantification", "Workflow orchestration"]
            },
            "gwas": {
                "template": "gwas_template.py",
                "description": "Genome-wide association studies",
                "features": ["Association testing", "Statistical analysis", "Visualization"]
            },
            "protein": {
                "template": "protein_template.py",
                "description": "Protein sequence analysis",
                "features": ["Structure analysis", "Sequence alignment", "Property calculation"]
            },
            "ml": {
                "template": "ml_template.py",
                "description": "Machine learning analysis",
                "features": ["Model training", "Feature selection", "Performance evaluation"]
            }
        }

    def generate_example(self, domain: str, name: str, description: str | None = None,
                        features: List[str] | None = None, template: str | None = None,
                        dry_run: bool = False) -> Dict[str, Any]:
        """Generate a new example file."""
        if domain not in self.domain_configs:
            available_domains = list(self.domain_configs.keys())
            raise ValueError(f"Unknown domain '{domain}'. Available: {available_domains}")

        # Get domain config
        config = self.domain_configs[domain]

        # Set defaults
        description = description or config["description"]
        features = features or config["features"]
        template = template or config["template"]

        # Template path
        template_path = self.templates_dir / template
        if not template_path.exists():
            raise FileNotFoundError(f"Template not found: {template_path}")

        # Output paths
        example_dir = self.examples_dir / domain
        example_file = example_dir / f"example_{name}.py"

        # Check if example already exists
        if example_file.exists():
            raise FileExistsError(f"Example already exists: {example_file}")

        # Load template
        with open(template_path, 'r', encoding='utf-8') as f:
            template_content = f.read()

        # Prepare template variables
        template_vars = {
            "domain": domain,
            "name": name,
            "description": description,
            "features": features,
            "analysis_type": self._get_analysis_type(domain, name),
            "imports": self._get_domain_imports(domain),
            "custom_imports": self._get_custom_imports(domain),
            "code_block": self._get_code_block(domain, name)
        }

        # Render template
        rendered_content = self._render_template(template_content, template_vars)

        if dry_run:
            return {
                "action": "dry_run",
                "example_file": str(example_file),
                "content": rendered_content,
                "template_vars": template_vars
            }

        # Create directory if needed
        example_dir.mkdir(parents=True, exist_ok=True)

        # Write example file
        with open(example_file, 'w', encoding='utf-8') as f:
            f.write(rendered_content)

        # Update README if it exists
        readme_path = example_dir / "README.md"
        if readme_path.exists():
            self._update_readme(readme_path, name, description)

        return {
            "action": "created",
            "example_file": str(example_file),
            "readme_updated": readme_path.exists(),
            "template_used": template
        }

    def _render_template(self, template: str, variables: Dict[str, Any]) -> str:
        """Render template with variable substitution."""
        # Simple template rendering (could be enhanced with Jinja2)
        result = template

        # Replace {{variable}} patterns
        for key, value in variables.items():
            if isinstance(value, list):
                # Handle list variables
                if key == "features":
                    # Special handling for features list
                    features_text = "\n".join(f"- {feature}" for feature in value)
                    result = result.replace("{% for feature in features -%}\n- {{ feature }}\n{% endfor %}", features_text)
                elif key == "imports":
                    imports_text = "\n".join(value)
                    result = result.replace("{% for import_line in imports -%}\n{{ import_line }}\n{% endfor %}", imports_text)
                else:
                    result = result.replace(f"{{{{ {key} }}}}", str(value))
            else:
                result = result.replace(f"{{{{ {key} }}}}", str(value))

        # Clean up any remaining template tags
        result = re.sub(r'{%.*?%}', '', result)
        result = re.sub(r'{{.*?}}', '', result)

        return result

    def _get_analysis_type(self, domain: str, name: str) -> str:
        """Get analysis type based on domain and name."""
        analysis_types = {
            "dna": {
                "sequences": "Sequence",
                "alignment": "Alignment",
                "phylogeny": "Phylogenetic",
                "population": "Population Genetics"
            },
            "rna": {
                "amalgkit": "Amalgkit Workflow",
                "quantification": "Expression Quantification"
            },
            "gwas": {
                "association": "Association Testing",
                "visualization": "GWAS Visualization"
            }
        }

        return analysis_types.get(domain, {}).get(name, "Analysis")

    def _get_domain_imports(self, domain: str) -> List[str]:
        """Get domain-specific imports."""
        imports_map = {
            "dna": ["from metainformant.dna import sequences, alignment"],
            "rna": ["from metainformant.rna import workflow"],
            "gwas": ["from metainformant.gwas import association"],
            "protein": ["from metainformant.protein import sequences"],
            "ml": ["from metainformant.ml import classification"],
            "core": []
        }
        return imports_map.get(domain, [])

    def _get_custom_imports(self, domain: str) -> str:
        """Get custom imports section."""
        return "# Additional imports as needed"

    def _get_code_block(self, domain: str, name: str) -> str:
        """Get the main code block for the example."""
        code_blocks = {
            "dna": '''
        # Perform sequence analysis
        if "{{analysis_type}}" == "Sequence":
            results_data[seq_id].update({
                "reverse_complement": sequences.reverse_complement(sequence),
                "gc_skew": sequences.gc_content(sequence) - 0.5,  # Simplified
                "complexity": len(set(sequence)) / len(sequence)  # Simplified
            })
        elif "{{analysis_type}}" == "Alignment":
            # Simple alignment example
            results_data[seq_id]["aligned_with_seq1"] = "Example alignment result"
        ''',
            "core": '''
        # Example core functionality
        results_data = {
            "config_loaded": True,
            "io_operations": ["load_json", "dump_json"],
            "logging_setup": True,
            "paths_validated": True
        }
        '''
        }

        return code_blocks.get(domain, '''
        # Add your example code here
        results_data = {
            "placeholder": "Replace with actual analysis results"
        }
        ''')

    def _update_readme(self, readme_path: Path, name: str, description: str) -> None:
        """Update domain README to include new example."""
        try:
            with open(readme_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # Add to examples list if there's a section
            if "## Examples" in content:
                example_entry = f"- `example_{name}.py`: {description}\n"
                content = content.replace("## Examples", f"## Examples\n{example_entry}", 1)

                with open(readme_path, 'w', encoding='utf-8') as f:
                    f.write(content)

        except Exception:
            # Don't fail if README update fails
            pass


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Generate METAINFORMANT examples")
    parser.add_argument("domain", help="Target domain (dna, rna, gwas, etc.)")
    parser.add_argument("name", help="Example name")
    parser.add_argument("--description", help="Brief description")
    parser.add_argument("--features", help="Comma-separated list of features")
    parser.add_argument("--template", help="Specific template to use")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be created")

    args = parser.parse_args()

    # Parse features
    features = None
    if args.features:
        features = [f.strip() for f in args.features.split(",")]

    # Create generator
    generator = ExampleGenerator()

    try:
        result = generator.generate_example(
            domain=args.domain,
            name=args.name,
            description=args.description,
            features=features,
            template=args.template,
            dry_run=args.dry_run
        )

        if args.dry_run:
            print("DRY RUN - Would create:")
            print(f"File: {result['example_file']}")
            print("\nContent preview:")
            print("-" * 40)
            print(result['content'][:500] + "..." if len(result['content']) > 500 else result['content'])
        else:
            print(f"✅ Created example: {result['example_file']}")
            if result.get('readme_updated'):
                print("✅ Updated README.md")
            print(f"Template used: {result['template_used']}")

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
