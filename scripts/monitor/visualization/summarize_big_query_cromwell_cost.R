#!/usr/bin/env Rscript
# "Script to pretty-display cost of runnning a cromwell job retrieved from Big Query"

################################################################################
library("optparse", warn.conflicts = F, quietly = T)

parser <- OptionParser()
parser <- add_option(parser, c("--wid"), type='character',
                     help="workflow id (root level) to query for cost")
parser <- add_option(parser, c("-b","--startDate"), type='character',
                     help="Year-Month-Day formatted date when the job started (allow some slack)")
parser <- add_option(parser, c("-e", "--endDate"), type='character',
                     help="Year-Month-Day formatted date when the job ended (allow some slack)")

parser <- add_option(parser, c("--plot"), type='character',
                     help="File name to store the cost plot")
parser <- add_option(parser, c("--width"), type='double',
                     default = 10,
                     help="Width of the plot")
parser <- add_option(parser, c("--height"), type='double',
                     default = 10,
                     help="Height of the plot")

parser <- add_option(parser, c("--json"), type='character',
                     default = NA,
                     help="File name to store the BQ result JSON file (optional)")
parser <- add_option(parser, c("--md"), type='character',
                     default = NA,
                     help="File name to store the detailed cost, formated in markdown (optional)")

################################################################################
summarize_cost <- function(workflow.id, start.date, end.date,
                           plot.pdf, width, height,
                           markdown.table.md = NA) {

    sql = get_query_sql(workflow.id, start.date, end.date)
    cost.df = run_query(sql)
    visualize(cost.df, plot.pdf, width, height, markdown.table.md)
}

################################################################################
# generate query SQL
get_query_sql <- function(workflow.id, start.date, end.date) {

    time.restriction = sprintf("_PARTITIONTIME BETWEEN TIMESTAMP(\'%s\') AND TIMESTAMP(\'%s\')",
                               start.date, end.date)

    wid.restriction = sprintf("label.value LIKE \"%%%s%%\"",
                              workflow.id)

    query = sprintf(
        "
        SELECT
          (SELECT value FROM UNNEST(labels) WHERE key = 'cromwell-workflow-id') AS workflow_id,
          (SELECT value FROM UNNEST(labels) WHERE key = 'cromwell-workflow-name') AS workflow_name,
          (SELECT value FROM UNNEST(labels) WHERE key = 'cromwell-sub-workflow-name') AS sub_workflow_id,
          (SELECT value FROM UNNEST(labels) WHERE key = 'wdl-task-name') AS task_name,
          (SELECT value FROM UNNEST(labels) WHERE key = 'wdl-call-alias') AS task_alias,
          (SELECT value FROM UNNEST(labels) WHERE key = 'goog-gke-node') AS node,

          sku.description AS sku_description,
          cost,
          cost_type,

          usage_start_time,
          usage_end_time,
          (SELECT value FROM UNNEST(system_labels) WHERE key = 'compute.googleapis.com/cores') AS cores,
          (SELECT value FROM UNNEST(system_labels) WHERE key = 'compute.googleapis.com/memory') AS memory,
          usage

        FROM
          `broad-dsde-methods.Methods_billing_dump.gcp_billing_export_v1_009C7D_923007_219A6F`

        LEFT JOIN
        UNNEST(labels) AS label

        WHERE
          %s
          AND cost > 0.0
          AND label.key IN (\"cromwell-workflow-id\",
                            \"cromwell-workflow-name\",
                            \"cromwell-sub-workflow-name\",
                            \"wdl-task-name\",
                            \"wdl-call-alias\",
                            \"goog-gke-node\")
          AND %s
        ",
        time.restriction,
        wid.restriction
    )
    query
}

################################################################################
# run the query, and slighly preprocess the returned data (data.frame)
run_query <- function(sql.query,
                      project.id = "broad-dsde-methods") {

    library("bigrquery", warn.conflicts = F, quietly = T)

    # TODO: this download approach is appropriate only for small tables
    # use API call to directly parse the result
    my.bq.table = bq_project_query(project.id, sql.query)
    df = bq_table_download(my.bq.table)

    df$"sku_description" = gsub(" running in Americas", "", df$"sku_description")

    df$"cost" = as.numeric(df$"cost")
    df = df[order(df$"cost", decreasing = T), ]

    df

    ## old solution that relies on forking out to a system command
    # temp.sql = tempfile(pattern = "bq.query", tmpdir = tempdir(), fileext = "sql")
    # if (is.null(query.json.file) || is.na(query.json.file)) {
    #     store.json = tempfile(pattern = "bq.result", tmpdir = tempdir(), fileext = "json")
    # } else {
    #     store.json = query.json.file
    # }


    # store.json
    # cmd = sprintf("bq --format=prettyjson query --nouse_legacy_sql \\
    #                   --flagfile=%s \\
    #                   > %s",
    #               temp.sql,
    #               store.json)
    # system(cmd)
}

################################################################################
visualize <- function(cost.df,
                      plot.pdf,
                      width, height,
                      markdown.table.md = NA) {

    library("knitr", warn.conflicts = F, quietly = T)
    library("dplyr", warn.conflicts = F, quietly = T)
    library("ggplot2", warn.conflicts = F, quietly = T)

    total.cost = sum(cost.df$"cost")

    # cut down, but other info maybe useful later
    presentation = cost.df[,c('cost', 'sku_description')]
    names.arr = apply(cost.df, 1,
                      function(e) {
                          if (is.na(e[['task_alias']])) e[['task_name']] else e[['task_alias']]
                      })
    presentation["task"] = names.arr
    names(presentation)= c("unitcost", "sku", "task")
    presentation$"sku" = as.factor(presentation$"sku")
    presentation$"task" = as.factor(presentation$"task")

    ##########
    # save markdown table, if requested
    if (!is.na(markdown.table.md) & !is.null(markdown.table.md)) {
        sink(markdown.table.md)
        md.table = kable(presentation, caption = sprintf("Total cost: %.2f", total.cost),
                         format = "markdown", row.names = F)
        print( md.table )
        sink()
    }

    ##########
    # save figures
    per.task <- presentation %>%
        group_by(task) %>%
        summarise(cost = sum(unitcost)) %>%
        mutate(share=cost/sum(cost)*100.0) %>%
        arrange(desc(cost))
    # # bar chart giving quick glance into which task costs how much percentage
    # # https://tinyurl.com/ya58xjtp
    # p = ggplot(per.task, aes("", share, fill = task)) +
    #     geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    #     coord_polar("y") +
    #     geom_text(aes(label = paste0(round(share), "%")),
    #               position = position_stack(vjust = 0.5)) +
    #     labs(x = NULL, y = NULL, fill = NULL,
    #          title = sprintf("Cost share per task class out of total cost: %.2f", total.cost)) +
    #     guides(fill = guide_legend(reverse = TRUE)) +
    #     theme_classic() +
    #     theme(axis.line = element_blank(),
    #           axis.text = element_blank(),
    #           axis.ticks = element_blank(),
    #           plot.title = element_text(hjust = 0.5))

    # real output, break down per task, then within each task, break down per SKU
    per.sku.per.task = NULL
    for (t in levels(presentation$"task")) {
        tmp <- presentation %>%
            filter(task == t) %>%
            group_by(sku) %>%
            summarise(cost = sum(unitcost)) %>%
            mutate(share=cost/sum(cost)*100.0) %>%
            arrange(desc(share))
        df = as.data.frame(tmp)
        df = subset(df, select=-c(cost))
        df$"sku" = as.factor(df$"sku")
        df$"task" = t
        if (is.null(per.sku.per.task)) {
            per.sku.per.task = df
        } else {
            per.sku.per.task = rbind.data.frame(per.sku.per.task, df)
        }
    }
    per.sku.per.task$"sku" = as.factor(per.sku.per.task$"sku")
    per.sku.per.task$"task" = factor(per.sku.per.task$"task",
                                     levels = as.character(per.task$task))

    cost.labeller = vector()
    for(t in per.task$"task") {
        s = sprintf("%s: $%.2f", t, per.task$"cost"[ per.task$"task" == t ])
        cost.labeller <- c(cost.labeller, s)
    }
    names(cost.labeller) = per.task$"task"
    n.row = floor( sqrt(length(levels(per.sku.per.task$"task"))) )
    q = ggplot(per.sku.per.task, aes("", share, fill = sku)) +
        geom_bar(color = "white", stat = "identity") +
        coord_polar("y") +
        geom_text(aes(label = paste0(round(share), "%")),
                  position = position_stack(vjust = 0.5)) +
        labs(x = NULL, y = NULL, fill = NULL,
             title = sprintf("Cost per task per SKU (tot. $%.2f)", total.cost)) +
        guides(fill = guide_legend(reverse = TRUE, ncol = 2)) +
        facet_wrap(facets=. ~ task, nrow = n.row,
                   labeller = as_labeller(cost.labeller))  +
        # theme_classic() +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position="bottom")
    # if the PDF looks cramed, strange, change the following parameters
    ggsave(plot.pdf, plot = q, width = width, height = height)
}
################################################################################

if (0 == length(commandArgs(trailingOnly=TRUE))) {
    parsed.args = parse_args(parser, args = c("--help"))
} else {
    parsed.args = parse_args(parser)
    summarize_cost(parsed.args$"wid", parsed.args$"startDate", parsed.args$"endDate",
                   parsed.args$"plot", parsed.args$"width", parsed.args$"height",
                   parsed.args$"md")
}
