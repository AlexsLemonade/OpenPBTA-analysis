# This is a Shiny web application designed to produce a d2b sunburst multilayer
# pie chart, which represents the histologies across samples within the given
# dataset.
#
# Chante Bethell for CCDL 2019

# Load/install packages.
if (!require("shinydashboard")) {
    install.packages("shinydashboard")
}
library(shinydashboard)
library(shiny)

# This will be removed once a Dockerfile is established
source("00-install-packages.R")

# Define UI for application that produces a d2b sunburst plot
ui <- dashboardPage(
    dashboardHeader(title = "Cancer Histologies"),
    dashboardSidebar(sidebarMenu(
        menuItem("Sunburst Plot", tabName = "sunbrstPlot")
    )),
    dashboardBody(
        tabBox(
            id = "sunbrstPlot",
            width = "20%",
            height = "1000px",
            sund2bOutput("Plot", height = "750", width = "100%")
        )
    )
)
# Define server logic required to produce the sunburst plot
server <- function(input, output) {
    df2 <- data.frame(readr::read_tsv(
        file.path(
            "..",
            "..",
            "data",
            "2019-07-26-pbta-histologies - 2019-07-26-including-germline.tsv"
        )
    ))
    
    sun_df <- df2 %>%
        dplyr::select(broad_histology, short_histology, disease_type_new) %>%
        dplyr::filter(!is.na(disease_type_new)) %>%
        dplyr::filter(!is.na(broad_histology)) %>%
        dplyr::filter(!is.na(short_histology)) %>%
        dplyr::group_by(broad_histology)
    
    sun_df_sub <- sun_df[, c("broad_histology", "short_histology",
                             "disease_type_new")]
    sun_df_sub$nodes <-
        paste0(sun_df_sub[[1]], ",", sun_df_sub[[2]], ",", sun_df_sub[[3]])
    
    # Make node patterns
    names(sun_df_sub) <- c("level1", "level2", "level3", "nodes")
    
    # Sum each unique node combination
    counts <- sun_df_sub %>%
        dplyr::group_by(nodes) %>%
        dplyr::summarise(size = n())
    
    counts <- counts[!grepl("^NA-", counts$nodes),]
    
    final <- merge(sun_df_sub, counts, all.x = T)
    col_order <- c("level1", "level2", "level3",
                   "size", "nodes")
    
    final <- final[, col_order]
    final$nodes <- NULL
    
    #colorblind-friendly color vector
    color <-
        c(
            "#386db0",
            "#E56E60",
            "#41A36D",
            "#F3E500",
            "#319106",
            "#0d76ff",
            "#C2B700",
            "#aac4e4",
            "#86cea6",
            "#999999",
            "#FFFBC9"
        )
    
    tm <-
        treemap::treemap(
            final,
            index = c("level1", "level2", "level3"),
            vSize = "size",
            vColor = color,
            draw = TRUE
        )
    
    tmnest <-
        d3r::d3_nest(tm$tm[, c("level1", "level2", "level3", "vSize")],
                     value_cols = c("vSize"))
    
    output$Plot <- renderSund2b({
        # Generate d2b sunburst plot using sund2b
        p <-
            sunburstR::sund2b(tmnest, colors = color, valueField = "vSize")
        print(p)
    })
}

# Run the application
shinyApp(ui = ui, server = server)
