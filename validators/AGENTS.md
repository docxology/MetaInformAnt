# Validators Module Technical Documentation

This document provides comprehensive technical documentation for the planned validator classes and interfaces in the MetaInformAnt validators module.

## Planned Components

### SchemaValidator Class

#### Constructor

```python
def __init__(self, schema_registry: Optional[Dict[str, Dict]] = None) -> None
```

**Parameters:**
- `schema_registry`: Optional initial schema registry dictionary

**Attributes:**
- `schemas`: Dictionary mapping schema names to schema definitions
- `cache`: LRU cache for compiled validation functions
- `strict_mode`: Boolean flag for strict validation behavior

#### validate

```python
def validate(self, data: Any, schema_name: str, **kwargs) -> ValidationResult
```

**Parameters:**
- `data`: Data to validate
- `schema_name`: Name of schema to validate against
- `**kwargs`: Additional validation options (strict_mode, fail_fast, etc.)

**Returns:**
- `ValidationResult`: Validation result with success status and error details

**Validation Process:**
1. Retrieve schema from registry
2. Compile schema if not cached
3. Validate data against schema
4. Collect validation errors
5. Return structured result

#### add_schema

```python
def add_schema(self, name: str, schema: Dict[str, Any]) -> None
```

**Parameters:**
- `name`: Unique schema identifier
- `schema`: JSON schema definition

**Effects:**
- Adds schema to registry
- Clears related cache entries
- Validates schema structure (if in strict mode)

#### remove_schema

```python
def remove_schema(self, name: str) -> bool
```

**Parameters:**
- `name`: Schema name to remove

**Returns:**
- `bool`: True if schema was removed, False if not found

#### list_schemas

```python
def list_schemas(self) -> List[str]
```

**Returns:**
- `List[str]`: List of registered schema names

#### get_schema

```python
def get_schema(self, name: str) -> Optional[Dict[str, Any]]
```

**Parameters:**
- `name`: Schema name to retrieve

**Returns:**
- `Optional[Dict[str, Any]]`: Schema definition or None if not found

### IntegrityValidator Class

#### Constructor

```python
def __init__(self, config: Optional[Dict[str, Any]] = None) -> None
```

**Configuration Parameters:**
- `algorithms`: List of supported hash algorithms (default: ['sha256', 'md5'])
- `encoding`: Data encoding for hashing (default: 'utf-8')
- `salt_length`: Salt length for salted hashing (default: 16)

#### validate_checksum

```python
def validate_checksum(self, data: Any, expected_hash: str, algorithm: str = 'sha256') -> IntegrityResult
```

**Parameters:**
- `data`: Data to validate
- `expected_hash`: Expected hash value
- `algorithm`: Hash algorithm to use

**Returns:**
- `IntegrityResult`: Result with validity, computed hash, and expected hash

#### generate_integrity_hash

```python
def generate_integrity_hash(self, data: Any, algorithm: str = 'sha256', salt: Optional[str] = None) -> str
```

**Parameters:**
- `data`: Data to hash
- `algorithm`: Hash algorithm ('sha256', 'md5', 'sha1', etc.)
- `salt`: Optional salt for salted hashing

**Returns:**
- `str`: Hexadecimal hash string

#### validate_file_integrity

```python
def validate_file_integrity(self, file_path: Union[str, Path], expected_hash: str, algorithm: str = 'sha256') -> IntegrityResult
```

**Parameters:**
- `file_path`: Path to file to validate
- `expected_hash`: Expected hash value
- `algorithm`: Hash algorithm to use

**Returns:**
- `IntegrityResult`: File integrity validation result

#### generate_salted_hash

```python
def generate_salted_hash(self, data: Any, algorithm: str = 'sha256') -> Tuple[str, str]
```

**Parameters:**
- `data`: Data to hash

**Returns:**
- `Tuple[str, str]`: (salt, salted_hash) pair

#### validate_salted_hash

```python
def validate_salted_hash(self, data: Any, salt: str, expected_hash: str, algorithm: str = 'sha256') -> bool
```

**Parameters:**
- `data`: Data to validate
- `salt`: Salt used during hashing
- `expected_hash`: Expected salted hash
- `algorithm`: Hash algorithm used

**Returns:**
- `bool`: True if hash matches

### ConsistencyValidator Class

#### Constructor

```python
def __init__(self, config: Optional[Dict[str, Any]] = None) -> None
```

**Configuration Parameters:**
- `max_datasets`: Maximum number of datasets for consistency checks
- `timeout`: Timeout for consistency operations (seconds)
- `cache_results`: Whether to cache consistency check results

#### validate_relationships

```python
def validate_relationships(self, datasets: List[Dict[str, Any]], relationships: List[Dict[str, Any]]) -> ConsistencyResult
```

**Parameters:**
- `datasets`: List of datasets to check
- `relationships`: List of relationship definitions

**Returns:**
- `ConsistencyResult`: Consistency check results

**Relationship Definition Format:**
```python
{
    'type': 'foreign_key',  # 'foreign_key', 'unique', 'referential'
    'primary_dataset': 'dataset1',
    'primary_key': 'id',
    'foreign_dataset': 'dataset2',
    'foreign_key': 'parent_id'
}
```

#### check_referential_integrity

```python
def check_referential_integrity(self, primary_data: Dict[str, Any], foreign_data: Dict[str, Any], key_mapping: Dict[str, str]) -> ConsistencyResult
```

**Parameters:**
- `primary_data`: Primary dataset (dict of records)
- `foreign_data`: Foreign dataset (dict of records)
- `key_mapping`: Mapping of primary keys to foreign keys

**Returns:**
- `ConsistencyResult`: Referential integrity check results

#### validate_unique_constraints

```python
def validate_unique_constraints(self, data: Dict[str, Any], unique_fields: List[str]) -> ConsistencyResult
```

**Parameters:**
- `data`: Dataset to check
- `unique_fields`: List of fields that should be unique

**Returns:**
- `ConsistencyResult`: Uniqueness validation results

#### validate_data_ranges

```python
def validate_data_ranges(self, data: Dict[str, Any], range_rules: Dict[str, Dict[str, Any]]) -> ConsistencyResult
```

**Parameters:**
- `data`: Dataset to validate
- `range_rules`: Range validation rules per field

**Range Rule Format:**
```python
{
    'field_name': {
        'min': 0,
        'max': 100,
        'type': 'numeric'  # 'numeric', 'date', 'string'
    }
}
```

**Returns:**
- `ConsistencyResult`: Range validation results

### BusinessRuleValidator Class

#### Constructor

```python
def __init__(self, rules: Optional[Dict[str, Callable]] = None) -> None
```

**Parameters:**
- `rules`: Optional initial rule registry

#### add_rule

```python
def add_rule(self, name: str, rule_func: Callable[[Any], ValidationResult]) -> None
```

**Parameters:**
- `name`: Rule identifier
- `rule_func`: Function that takes data and returns ValidationResult

#### validate_business_rules

```python
def validate_business_rules(self, data: Any, rule_names: List[str]) -> List[ValidationResult]
```

**Parameters:**
- `data`: Data to validate
- `rule_names`: List of rule names to apply

**Returns:**
- `List[ValidationResult]`: Results for each applied rule

#### remove_rule

```python
def remove_rule(self, name: str) -> bool
```

**Parameters:**
- `name`: Rule name to remove

**Returns:**
- `bool`: True if rule was removed

## Data Structures

### ValidationResult

```python
@dataclass
class ValidationResult:
    valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
```

### IntegrityResult

```python
@dataclass
class IntegrityResult:
    valid: bool
    computed_hash: str
    expected_hash: str
    algorithm: str
    data_size: int
    validation_time: float
```

### ConsistencyResult

```python
@dataclass
class ConsistencyResult:
    consistent: bool
    violations: List[Dict[str, Any]] = field(default_factory=list)
    checked_items: int = 0
    validation_time: float = 0.0
```

## Exception Classes

### ValidationError

```python
class ValidationError(MetaInformAntError):
    def __init__(
        self,
        message: str,
        validator_type: Optional[str] = None,
        schema_name: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None
    ) -> None
```

**Additional Attributes:**
- `validator_type`: Type of validator that failed
- `schema_name`: Schema name if applicable

### IntegrityError

```python
class IntegrityError(MetaInformAntError):
    def __init__(
        self,
        message: str,
        algorithm: Optional[str] = None,
        expected_hash: Optional[str] = None,
        computed_hash: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None
    ) -> None
```

**Additional Attributes:**
- `algorithm`: Hash algorithm used
- `expected_hash`: Expected hash value
- `computed_hash`: Computed hash value

### ConsistencyError

```python
class ConsistencyError(MetaInformAntError):
    def __init__(
        self,
        message: str,
        violation_count: int = 0,
        datasets_checked: Optional[List[str]] = None,
        context: Optional[Dict[str, Any]] = None
    ) -> None
```

**Additional Attributes:**
- `violation_count`: Number of consistency violations found
- `datasets_checked`: List of dataset names checked

## Configuration Schema

```python
ValidatorConfig = Dict[str, Any]  # Validator configuration
SchemaDefinition = Dict[str, Any]  # JSON schema definition
ValidationRule = Callable[[Any], ValidationResult]  # Validation rule function
RelationshipDefinition = Dict[str, Any]  # Relationship definition for consistency checks
RangeRule = Dict[str, Dict[str, Any]]  # Range validation rules
```

## Performance Characteristics

### Validation Performance

- **Schema Validation**: O(n) where n is data size, with caching for repeated schemas
- **Integrity Checks**: O(m) where m is data size, dominated by hash computation
- **Consistency Checks**: O(d × r) where d is dataset count, r is relationship count
- **Business Rules**: Depends on rule complexity, typically O(n) for data size n

### Optimization Strategies

- **Schema Caching**: Compiled schemas cached with LRU eviction
- **Incremental Validation**: Validate only changed portions of data
- **Parallel Processing**: Multi-threaded validation for large datasets
- **Batch Operations**: Process multiple validations together

## Integration Patterns

### Agent Integration

```python
class ValidationAgent(BaseAgent):
    def get_capabilities(self) -> List[str]:
        return ['validate_schema', 'check_integrity', 'verify_consistency']

    def execute(self, task: Task) -> TaskResult:
        if task.operation == 'validate_schema':
            return self._validate_schema_task(task)
        elif task.operation == 'check_integrity':
            return self._check_integrity_task(task)
        # ... other validation operations
```

### Orchestrator Integration

```python
# Register validation capabilities
orchestrator.register_agent(validation_agent, [
    'validate_schema',
    'check_integrity',
    'verify_consistency',
    'validate_business_rules'
])

# Coordinate validation tasks
validation_task = Task(
    operation='validate_schema',
    parameters={
        'data': dataset,
        'schema_name': 'user_data_schema'
    }
)
result = orchestrator.coordinate_task(validation_task)
```

This module will be implemented following MetaInformAnt's modular architecture with comprehensive error handling, logging, and performance monitoring.

