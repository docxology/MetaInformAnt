# AI Agents Documentation

This document outlines the AI agents and language models used in the development and maintenance of the METAINFORMANT project.

## Project Overview

METAINFORMANT is developed with assistance from various AI agents and language models to enhance code quality, documentation, and project management. This collaborative approach leverages AI capabilities for:

- **Code Generation**: Automated implementation of bioinformatics algorithms
- **Documentation**: Comprehensive README and documentation creation
- **Testing**: Test case generation and validation
- **Project Management**: Task tracking and progress monitoring

## AI Agents Used

### Primary Development Agent
**Code Assistant Agent** - Cursor's AI coding assistant
- **Model**: grok-code-fast-1
- **Purpose**: Real-time code assistance, file editing, and project management
- **Capabilities**:
  - Code generation and refactoring
  - Documentation writing and validation
  - Test case generation and validation
  - Bug detection and fixing
  - Project structure optimization
  - Multi-omic bioinformatics algorithm implementation
  - Integration of scientific computing libraries

### Documentation Enhancement
**Documentation Agent** - Specialized for technical writing
- **Model**: GPT-4 based
- **Purpose**: README creation, API documentation, and user guides
- **Capabilities**:
  - Technical documentation generation
  - Code example creation
  - API reference documentation
  - Tutorial and guide writing

### Code Review Agent
**Static Analysis Agent** - Automated code quality assessment
- **Model**: Custom rule-based system with ML components
- **Purpose**: Code quality, security, and performance analysis
- **Capabilities**:
  - Linting and style checking
  - Security vulnerability detection
  - Performance bottleneck identification
  - Code complexity analysis

## AI Integration Workflow

### Development Process
1. **Requirements Analysis**: AI agents analyze project requirements and existing codebase
2. **Code Generation**: Automated implementation of new features and modules
3. **Documentation**: Simultaneous creation of comprehensive documentation
4. **Testing**: Automated test case generation and validation
5. **Review**: AI-assisted code review and quality assurance

### Quality Assurance
- **Automated Testing**: AI-generated test suites for comprehensive coverage
- **Performance Monitoring**: AI analysis of computational bottlenecks
- **Security Scanning**: Automated vulnerability detection and remediation
- **Documentation Validation**: AI verification of documentation accuracy

## AI-Generated Content

### Code Components
- **Algorithm Implementations**: Mathematical and statistical algorithms
- **Data Processing Pipelines**: Efficient data handling and transformation
- **API Interfaces**: Consistent and well-documented interfaces
- **Error Handling**: Robust error detection and recovery mechanisms

### Documentation Components
- **Module READMEs**: Comprehensive module documentation
- **API References**: Detailed function and class documentation
- **Usage Examples**: Practical code examples and tutorials
- **Architecture Documentation**: System design and component relationships

### Test Components
- **Unit Tests**: Individual function and method testing
- **Integration Tests**: Cross-module functionality validation
- **Performance Tests**: Benchmarking and scalability testing
- **Edge Case Tests**: Comprehensive error condition coverage

## Ethical Considerations

### Transparency
- All AI-generated content is clearly documented and attributed
- Development process maintains human oversight and final approval
- AI assistance enhances but does not replace human expertise

### Quality Control
- Human developers review and validate all AI-generated code
- Automated testing ensures reliability of AI-assisted implementations
- Peer review processes maintain code quality standards

### Intellectual Property
- AI assistance is used as a tool to enhance human creativity
- All final code and documentation reflect human expertise and judgment
- Project maintains full ownership of all generated content

## Best Practices

### AI Integration Guidelines
1. **Human Oversight**: All AI-generated content requires human review
2. **Transparency**: Clearly mark AI-assisted sections in documentation
3. **Validation**: Comprehensive testing of AI-generated code
4. **Ethical Use**: Responsible application of AI technologies

### Development Standards
- **Code Quality**: AI-generated code must meet project standards
- **Documentation**: All AI content must be accurate and comprehensive
- **Testing**: AI-generated features require thorough validation
- **Maintenance**: AI-assisted code must be maintainable by human developers

### Cursor Rules Compliance
- **Follow `.cursorrules`**: All AI agents must adhere to the main `.cursorrules` file
- **Module-Specific Rules**: Consult `cursorrules/<module>.cursorrules` for domain-specific patterns
- **See `cursorrules/README.md`**: For guidance on using modular cursorrules
- **Key Requirements**:
  - Write outputs to `output/` by default
  - Use `config/` with env overrides for configuration
  - No mocks in tests (real implementations only)
  - Use `metainformant.core` utilities for I/O, logging, paths
  - Update existing docs, never create root-level docs

## Future AI Integration

### Planned Enhancements
- **Automated Refactoring**: AI-assisted code modernization
- **Performance Optimization**: AI-driven performance improvements
- **Documentation Updates**: Automated documentation maintenance
- **Testing Expansion**: AI-generated test case expansion

### Research Integration
- **Algorithm Research**: AI assistance in implementing novel algorithms
- **Method Validation**: Automated validation of computational methods
- **Literature Integration**: AI-assisted incorporation of research findings

## Documentation

For detailed documentation about AI contributions to specific modules, organized by category:

### Repository-Level Documentation
- **Main Source Development**: [`src/metainformant/AGENTS.md`](src/metainformant/AGENTS.md) - Overall AI assistance in source code development
- **Source Organization**: [`src/AGENTS.md`](src/AGENTS.md) - Source code infrastructure and organization
- **Configuration Management**: [`config/AGENTS.md`](config/AGENTS.md) - Configuration system development
- **General Documentation**: [`docs/AGENTS.md`](docs/AGENTS.md) - AI assistance in documentation development

Note: The `output/` directory is ephemeral and contains generated analysis results. Documentation about output structure is in the `.cursorrules` file and module-specific documentation.

### Source Module Documentation (Implementation)
- **Core Utilities**: [`src/metainformant/core/AGENTS.md`](src/metainformant/core/AGENTS.md) - Core infrastructure and shared utilities
- **DNA Analysis**: [`src/metainformant/dna/AGENTS.md`](src/metainformant/dna/AGENTS.md) - DNA sequence analysis and genomics
- **RNA Analysis**: [`src/metainformant/rna/AGENTS.md`](src/metainformant/rna/AGENTS.md) - RNA transcriptomic analysis and workflow orchestration
- **RNA Workflow Steps**: [`src/metainformant/rna/steps/AGENTS.md`](src/metainformant/rna/steps/AGENTS.md) - Modular RNA workflow step implementations
- **GWAS Module**: [`src/metainformant/gwas/AGENTS.md`](src/metainformant/gwas/AGENTS.md) - Genome-wide association studies implementation
- **Protein Analysis**: [`src/metainformant/protein/AGENTS.md`](src/metainformant/protein/AGENTS.md) - Protein sequence and structure analysis
- **Mathematical Biology**: [`src/metainformant/math/AGENTS.md`](src/metainformant/math/AGENTS.md) - Mathematical and theoretical biology
- **Machine Learning**: [`src/metainformant/ml/AGENTS.md`](src/metainformant/ml/AGENTS.md) - Machine learning for biological data
- **Network Analysis**: [`src/metainformant/networks/AGENTS.md`](src/metainformant/networks/AGENTS.md) - Biological network analysis
- **Multi-Omics**: [`src/metainformant/multiomics/AGENTS.md`](src/metainformant/multiomics/AGENTS.md) - Multi-omic data integration
- **Single-Cell Genomics**: [`src/metainformant/singlecell/AGENTS.md`](src/metainformant/singlecell/AGENTS.md) - Single-cell RNA sequencing analysis
- **Quality Control**: [`src/metainformant/quality/AGENTS.md`](src/metainformant/quality/AGENTS.md) - Data quality assessment
- **Visualization**: [`src/metainformant/visualization/AGENTS.md`](src/metainformant/visualization/AGENTS.md) - Plotting and visualization utilities
- **Simulation**: [`src/metainformant/simulation/AGENTS.md`](src/metainformant/simulation/AGENTS.md) - Synthetic data generation and modeling
- **Ontology**: [`src/metainformant/ontology/AGENTS.md`](src/metainformant/ontology/AGENTS.md) - Functional annotation and ontologies
- **Phenotype**: [`src/metainformant/phenotype/AGENTS.md`](src/metainformant/phenotype/AGENTS.md) - Phenotypic trait analysis
- **Epigenome**: [`src/metainformant/epigenome/AGENTS.md`](src/metainformant/epigenome/AGENTS.md) - Epigenetic modification analysis
- **Ecology**: [`src/metainformant/ecology/AGENTS.md`](src/metainformant/ecology/AGENTS.md) - Ecological metadata and community analysis
- **Information Theory**: [`src/metainformant/information/AGENTS.md`](src/metainformant/information/AGENTS.md) - Information-theoretic analysis
- **Life Events**: [`src/metainformant/life_events/AGENTS.md`](src/metainformant/life_events/AGENTS.md) - Life course event analysis

### Documentation Module Files (User Documentation)
- **Core Documentation**: [`docs/core/AGENTS.md`](docs/core/AGENTS.md) - Core utilities documentation development
- **DNA Documentation**: [`docs/dna/AGENTS.md`](docs/dna/AGENTS.md) - DNA analysis documentation development
- **GWAS Documentation**: [`docs/gwas/AGENTS.md`](docs/gwas/AGENTS.md) - GWAS module documentation development
- **Mathematical Biology Documentation**: [`docs/math/AGENTS.md`](docs/math/AGENTS.md) - Mathematical biology documentation development
- **Machine Learning Documentation**: [`docs/ml/AGENTS.md`](docs/ml/AGENTS.md) - Machine learning documentation development
- **Network Analysis Documentation**: [`docs/networks/AGENTS.md`](docs/networks/AGENTS.md) - Network analysis documentation development
- **Quality Control Documentation**: [`docs/quality/AGENTS.md`](docs/quality/AGENTS.md) - Quality control documentation development
- **Ontology Documentation**: [`docs/ontology/AGENTS.md`](docs/ontology/AGENTS.md) - Ontology documentation development
- **Protein Documentation**: [`docs/protein/AGENTS.md`](docs/protein/AGENTS.md) - Protein analysis documentation development
- **Phenotype Documentation**: [`docs/phenotype/AGENTS.md`](docs/phenotype/AGENTS.md) - Phenotype documentation development
- **Epigenome Documentation**: [`docs/epigenome/AGENTS.md`](docs/epigenome/AGENTS.md) - Epigenome documentation development
- **Ecology Documentation**: [`docs/ecology/AGENTS.md`](docs/ecology/AGENTS.md) - Ecology documentation development
- **Visualization Documentation**: [`docs/visualization/AGENTS.md`](docs/visualization/AGENTS.md) - Visualization documentation development
- **Single-Cell Documentation**: [`docs/singlecell/AGENTS.md`](docs/singlecell/AGENTS.md) - Single-cell documentation development
- **Simulation Documentation**: [`docs/simulation/agents.md`](docs/simulation/agents.md) - Simulation documentation development
- **Multi-Omics Documentation**: [`docs/multiomics/AGENTS.md`](docs/multiomics/AGENTS.md) - Multi-omics documentation development
- **RNA Documentation**: [`docs/rna/AGENTS.md`](docs/rna/AGENTS.md) - RNA analysis documentation development
- **Amalgkit Documentation**: [`docs/rna/amalgkit/AGENTS.md`](docs/rna/amalgkit/AGENTS.md) - Amalgkit integration documentation
- **Amalgkit Steps Documentation**: [`docs/rna/amalgkit/steps/AGENTS.md`](docs/rna/amalgkit/steps/AGENTS.md) - Amalgkit workflow steps documentation
- **RNA Examples Documentation**: [`docs/rna/examples/AGENTS.md`](docs/rna/examples/AGENTS.md) - RNA workflow example guides

### Scripts and Testing
- **RNA Scripts**: [`scripts/rna/AGENTS.md`](scripts/rna/AGENTS.md) - RNA workflow scripts development
- **Amalgkit Scripts**: [`scripts/rna/amalgkit/AGENTS.md`](scripts/rna/amalgkit/AGENTS.md) - Amalgkit script development
- **Core Scripts**: [`scripts/core/AGENTS.md`](scripts/core/AGENTS.md) - Core utility scripts development
- **Test Data**: [`tests/data/AGENTS.md`](tests/data/AGENTS.md) - Test data organization documentation

## Contact and Support

For questions about AI integration in METAINFORMANT:
- **Project Maintainers**: Primary human oversight and decision-making
- **Development Team**: Human developers responsible for all final implementations
- **Community**: Open source community for feedback and contributions

---

*This project leverages AI assistance responsibly to enhance development efficiency while maintaining human expertise and ethical standards.*
