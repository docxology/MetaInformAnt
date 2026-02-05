#!/usr/bin/env python3
"""{{description}}

This example demonstrates:
{% for feature in features -%}
- {{ feature }}
{% endfor %}

Usage:
    python examples/{{domain}}/example_{{name}}.py

Expected output:
    output/examples/{{domain}}/{{name}}_results.json
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core import io

{% if imports %}
{% for import_line in imports -%}
{{ import_line }}
{% endfor %}
{% endif %}

{{custom_imports}}


def main():
    """Main function demonstrating {{description.lower()}}."""
    print("{{description}} Example")
    print("=" * 40)

    # Create output directory
    output_dir = Path("output/examples/{{domain}}")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "{{name}}_results.json"

    try:
        {{code_block}}

        # Save results
        results = {
            "example": "{{name}}",
            "domain": "{{domain}}",
            "description": "{{description}}",
            "results": results_data
        }

        io.dump_json(results, output_file)
        print(f"✅ Results saved to: {output_file}")

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
