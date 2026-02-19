# LIFE_EVENTS

## Overview
Life events and trajectory analysis module for METAINFORMANT.

## 📦 Contents
- **[analysis/](analysis/)**
- **[core/](core/)**
- **[models/](models/)**
- **[visualization/](visualization/)**
- **[workflow/](workflow/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

```mermaid
graph TD
    subgraph "Life Events Module"
        CE[core/events.py] --> |Event, EventSequence| DS[Data Structures]
        CE --> |EventDatabase| DB[Event Storage]
        CC[core/config.py] --> CF[Configuration]

        ME[models/embeddings.py] --> |learn_event_embeddings| EM[Word2Vec-style Embeddings]
        ME --> |biological_embedding| BE[Domain Embeddings]

        AI[analysis/interpretability.py] --> |attention_weights| AW[Model Interpretation]

        SV[survival/time_to_event.py] --> |kaplan_meier_estimator| KM[Kaplan-Meier]
        SV --> |cox_ph_model| COX[Cox PH]
        SV --> |competing_risks| CR[Competing Risks]

        WF[workflow/workflow.py] --> |analyze_life_course| LC[Life Course Analysis]
        WF --> |compare_populations| CP[Population Comparison]

        VZ[visualization/] --> VP[Event Plots]
    end
```

## Usage
Import module:
```python
from metainformant.life_events import ...
```
