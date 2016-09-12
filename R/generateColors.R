# Generate biome colors for plots

#' @export
generate.colors <- function(models){
  
  library(RColorBrewer)
  library(data.table)
  
  pft.biome <- try.pfts[,pft.num := sprintf("pft.%02.f",as.numeric(pft.factor))][,list(biome = unique(biome)),pft.num][,biome.num := as.numeric(as.factor(biome))][1:35];setnames(pft.biome,"biome","Biome")
  
  pal <- brewer.pal(length(unique(pft.biome[,Biome])),"Set1")
  
  pft.biome <- rbind(pft.biome, data.table(pft.num = "global", Biome = "global", biome.num = 0))
  
  color.dt <- data.table(Model = unique(models))[,pft.num :=  sapply(X = unique(models), FUN = function(x) tail(unlist(strsplit(x,"\\.na.")),1))]
  color.dt <- merge(color.dt,pft.biome,by = "pft.num", all.x = T)[,Color:= unlist(lapply(biome.num, function(x) ifelse(x == 0,"black", pal[x])))]
  
  return(color.dt[,.(Model,Biome,Color,pft.num)])
}
