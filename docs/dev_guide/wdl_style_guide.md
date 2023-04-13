# WDL Style Guide
This document describes the style guide for writing WDL workflows and tasks for the pipelines.
The guide will provide a list of best practices for writing WDL workflows, and will also 
provide a list of common mistakes to avoid. The guide is for those wanting to contribute
Dockerfiles to the pipelines repository.

## WDL Structure: 
This section describes the overall structure of a WDL workflow, including the input and output declarations, the task section, and the workflow section.

* A common WDL file contains both a workflow block and one or more task blocks. In the pipeline repository, 
it is recommended to place task blocks in a separate wdl file from the wdl file containing the workflow block.
This is to help keep the workflow WDLs more modular, readable, and maintainable.
* The WDL file containing the workflow block should include meta and parameter_meta 
sections. The meta section should include a description of the workflow, and the 
parameter_meta section should include a description of the workflow inputs. 
* For wdl files composed of task blocks it is _recommended_ to include meta and meta_parameters sections within the task blocks.

## Task Definitions: 
This section provides guidance on how to define tasks in WDL, including best practices for naming tasks, specifying inputs and outputs, and defining command scripts.

The order of the sections in a task block should be as follows:
1. meta
2. parameter_meta
3. input
4. command
5. output
6. runtime

The following is an example of a task block in WDL:
```WDL

task TaskName {
  meta {
    description: "A description of the task"
  }
  parameter_meta {
    input1: "A description of the input1"
    input2: "A description of the input2"
  }
  input {
    # task inputs
  }
  command {  
    # task command
  }
  output {
    # task outputs
  }
  runtime {
    # task runtime
  }
}
```

### Task Runtime:
The WDLs in the pipeline repository often use a runtime object from the struct.wdl in the task block. Using 
the struct isn't required but is recommended. 

The following is an example of a runtime struct used in a task block:
```WDL
import "structs.wdl" as structs
    
task TaskName {
    input {
        ...
        RuntimeAttr? runtime_attr_override
    }
    ...
    output {...}
        
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
```

## Workflow Definitions: 
This section covers best practices for defining workflows in WDL, including how to specify dependencies between tasks, how to handle errors and exceptions, and how to define scatter and gather blocks._

The order of the sections in a task block should be as follows:
1. meta
2. parameter_meta
3. input
4. call
5. output

The following is an example of a workflow block in WDL:
```WDL
workflow WorkflowName {
  meta {
      description: "A description of the workflow"
  }
  parameter_meta {
        input1: "A description of the input1"
        input2: "A description of the input2"
  }
  input {
    # workflow inputs
  }
  call task_name {
    # task inputs
  }

  output {
    # workflow outputs
  }
}
```

