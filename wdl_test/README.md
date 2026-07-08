# Pipeline Testing

## Testing Schedule 
- Run all available WDL workflows once a week
- Run a handful of important workflows when a push is made to feature branch


Note: Call caching is enabled so that if a wdl hasn't changed then its tests will not run. 

### Questions

- If an utils task is changed, theoretically, it affects all pipelines that use the utils task. Do we want to test all pipelines? Or do we want to have multiple tests for that task covering all scenarios? 
- This *should* be handled by the feature branch test. The test will attempt to run a workflow, if that workflow imports a util task that has been changed then cromwell should not call cache that imported workflow/task and instead and run the workflow from scratch.
