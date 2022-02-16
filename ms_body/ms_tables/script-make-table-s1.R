
\newpage

**Table S1**  
  
  \  

```{r table_S1, echo = FALSE}
# Import data to make Table S1
tableS1_data <- readr::read_csv(here::here("./data/data_clean/table2_data.csv")) %>%
  janitor::clean_names()


# Create Table S1
tableS1 <- tableS1_data %>%
  knitr::kable(
    format = "latex",
    align = "l",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    escape = FALSE,
    col.names = c(
      "Model",
      "RM",
      "FC",
      "AUC\\textsubscript{test}",
      "AUC\\textsubscript{diff}",
      "OR\\textsubscript{10}",
      "AIC\\textsubscript{c}"
    )
  ) %>%
  kableExtra::kable_styling(
    position = "left",
    latex_options = c("striped", "repeat_header")
    #,    stripe_color = "gray!15"
  )

# Display Table S1
tableS1
```



