library(ggplot2)
library(foreach)
library(Cairo)

for(i in 1){
  jpeg(paste0("~/Documents/BFX_proj/_small_tools/test/parallel_ggplot/", i, ".jpg"), height = 300, width = 300, type = "cairo", res = 150)
  print(
    ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
      geom_point() +
      ggtitle(i)
  )
  dev.off()
}


