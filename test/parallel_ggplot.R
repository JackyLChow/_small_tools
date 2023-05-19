library(ggplot2)
library(parallel)
library(foreach)


for(i in 1){
  jpeg(paste0("~/Documents/BFX_proj/_small_tools/test/parallel_ggplot/", i, ".jpg"), height = 300, width = 300, res = 150)
  ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
    geom_point() +
    ggtitle(i)
  dev.off()
}

system.time(
  for(i in 1:10){
    jpeg(paste0("~/Documents/BFX_proj/_small_tools/test/parallel_ggplot/", i, ".jpg"), height = 300, width = 300, res = 150)
    print(
      ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
        geom_point() +
        ggtitle(i) 
    )
    dev.off()
  }
)

plots <- list()
system.time(
  for(i in 1:1000){
    plots[[i]] <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
      geom_point() +
      ggtitle(i)
  }
)

system.time(
  plots <- mclapply(1:1000,
                    function(i){
                      ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
                        geom_point() +
                        ggtitle(i)
                    },
                    mc.cores = 10
  )
)

system.time(
  for(i in 1:length(plots)){
    jpeg(paste0("~/Documents/BFX_proj/_small_tools/test/parallel_ggplot/", i, ".jpg"), height = 300, width = 300, res = 150)
    print(
      plots[[i]]
    )
  }
)
