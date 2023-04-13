# WDL Style Guide
This document describes the style guide for writing WDL workflows and tasks for the pipelines.
The guide will provide a list of best practices for writing WDL workflows, and will also 
provide a list of common mistakes to avoid. The guide is for those wanting to contribute
Dockerfiles to the pipelines repository.

## Syntax and Formatting: 
_This section would cover the basic syntax and formatting rules for writing WDL workflows, including conventions for naming variables and tasks, indenting, and commenting._

## WDL Structure: 
_This section would describe the overall structure of a WDL workflow, including the input and output declarations, the task section, and the workflow section._

* A common WDL file contains both a workflow block and one or more task blocks. In the pipeline repository, 
it is recommended to place task blocks in a separate wdl file from the wdl workflow 
calling the task. This is to help keep the workflow wdls more modular, readable, and maintainable.
* The WDL file containing the workflow block should include meta and meta_parameters 
sections. The meta section should include a description of the workflow, and the 
meta_parameters section should include a description of the workflow inputs. 
* For wdl files composed of task blocks it is _recommended_ to include meta and meta_parameters sections within the task blocks.

### Workflow basic structure
```WDL
workflow workflow_name {
  meta {
      description: "A description of the workflow"
  }
  meta_parameters {
        input1: "A description of the input1"
        input2: "A description of the input2"
  }
  input {
    # workflow inputs
  }

  # task calls

  output {
    # workflow outputs
  }
}
```



## Task Definitions: 
_This section would provide guidance on how to define tasks in WDL, including best practices for naming tasks, specifying inputs and outputs, and defining command scripts._

## Workflow Definitions: 
_This section would cover best practices for defining workflows in WDL, including how to specify dependencies between tasks, how to handle errors and exceptions, and how to define scatter and gather blocks._

## Function Definitions: 
_This section would cover how to define reusable functions in WDL, including how to pass arguments and return values, and how to incorporate functions into workflows._

## Conventions and Best Practices: 
_This section would provide guidance on conventions and best practices for writing clear, maintainable, and portable WDL workflows, including recommendations for naming conventions, commenting, error handling, and testing._

Examples: This section would include examples of well-written WDL workflows, tasks, and functions, along with explanations of the rationale behind their design choices and the benefits they provide.