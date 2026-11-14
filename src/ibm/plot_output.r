#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}

find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}
# use tidyverse's sym() to transform strings to symbols 
transform.sym  <- function(x) 
{
    if (!is.null(x))
    {
        sym(x)
    }
}

jsonstuff <- '[
    {
        "xvar" : "time",
        "yvar" : "meanv"
    },
    {
        "xvar" : "time",
        "yvar" : "varv"
    },{
        "xvar" : "time",
        "yvar" : ["mean_baseline_influx","mean_stress_influx","mean_vigilance_influx"]
    },{
        "xvar" : "time",
        "yvar" : ["var_baseline_influx","var_stress_influx","var_vigilance_influx"]
    },{
        "xvar" : "time",
        "yvar" : "mean_removal"
    },{
        "xvar" : "time",
        "yvar" : ["mean_damage","mean_stress_hormone","mean_group_V"]
    }]
'
if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)

data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-1
        ,col_names=T)

# get the parameters
data.tibble.params <- read_delim(file=file.name
        ,delim=";"
        ,skip=param.line
        ,col_names=c("name","value")
        )

# transpose the tibble with the parameters
params <- data.tibble.params %>% pivot_wider(
        names_from = name
        ,values_from = value)

if (nrow(data.tibble.params) > 10000)
{
    data.tibble <- data.tibble %>% filter(generation %% 100 == 0)
}

plot.structure <- fromJSON(jsonstuff, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list()

plot.list.idx <- 1

for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        # get rid of aes_string like this: 
        # https://stackoverflow.com/questions/74414272/how-to-replace-the-deprecated-ggplot2-function-aes-string-accepting-an-arbitrar/74414389#74414389 

        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=yvar_values),
                FUN=transform.sym)

        # aes arguments for the geom_line() call
        line_args <- lapply(X=list(
                        colour=yvar_name),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=sub.data
                ,mapping=aes(!!!plot_args)) + 
                    geom_line(mapping=aes(!!!line_args))
    } else {
        
        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=plot.structure[[plot_struct_idx]]$yvar
                        ),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=data.tibble
                ,mapping=aes(!!!plot_args)) + geom_line()
    }

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
} # end for


# ok now add a file that simulates a stress response
tmax_iter <- 500
t_stressor <- 250

data_stress <- data.frame(
        t = seq(0,tmax_iter),
        stress_hormone = NA) %>% 
    mutate(
            stressor = if_else(t == t_stressor,1,0)
            )

print(params)

# start the stress hormone at its optimum
hormone_optimum <- as.numeric(params["theta_hormone"])
data_stress[1,"stress_hormone"] <- hormone_optimum

baseline_influx <- data.tibble %>% tail(n = 1) %>% pull(mean_baseline_influx)
stress_influx <- data.tibble %>% tail(n = 1) %>% pull(mean_stress_influx)
g <- as.numeric(params["g"])

for (t in 2:max(data_stress$t)) 
{
    hormone_t = data_stress[t-1,] %>% pull(stress_hormone)
    stressor = data_stress[t,] %>% pull(stressor)
    hormone_tplus1 = hormone_t * (1.0 - g) + 
        baseline_influx + 
        stress_influx * stressor

    data_stress[t,"stress_hormone"] = hormone_tplus1
}

p_dynamics <- ggplot(
        data = data_stress,
        mapping = aes(x = t, y = stress_hormone)) +
    geom_line()

plot.list[[plot.list.idx]] <- p_dynamics

title <- ""

if (exists("params") && "s_np" %in% names(params)) {
    title <- paste0(
	"npatches: ", params[["npatches"]], ", n: ", params[["n"]], "\n",
    "maxtime: ", params[["max_time"]], ", s_np: ", params[["s_np"]], "\n",
    "s_p: ", params[["s_p"]], ", m_d: ", params[["md"]], "\n",
    "m_v: ", params[["mv"]], ", p_mig: ", params[["p_mig"]], "\n",
    "p_attack: ", params[["p_attack"]], ", init_v: ", params[["init_v"]], "\n",
    "init_stress_hormone: ", params[["init_stress_hormone_level"]], "\n",
    "damage removal: ", params[["g"]], ", hormone damage incr: ", params[["k"]]
)

}

wrap_plots(plot.list,ncol=1) + 
	plot_annotation(title=title) +
	theme (plot.title = element_text(size=10))

file.name <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name,height= 3 * plot.structure.l)






