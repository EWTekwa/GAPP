library(spocc)


one <- occ(query = 'Enophrys bison', from = 'obis', limit = 5)
dd <- one$obis
dd[2]

#sometime in the 80s Clupea was split into different species! this presents a problem - the proble likely will persist across other genera

a1 <- occ(query = c('Enophrys bison','Clupea pallasii','Clupea harengus harengus', 'Rhacochilus vacca', 'Hippoglossus hippoglossus', 'Ammodytes japonicus', 'Ammodytes hexapterus'), 
          date = c('1940-08-01', '2020-08-31'),
          from = 'obis', limit = 5, geometry = geo$geometry)

b1 <- occ( 
          date = c('1945-08-01', '2020-08-31'),
          from = 'obis', limit = 10000, geometry = geo$geometry)

a2 <- a1$obis

a2 <- occ2df(b1)

print(a1$obis)



a1 <- occ(query = c('Enophrys bison','Clupea harengus', 'Rhacochilus vacca', 'Hippoglossus hippoglossus'), from = 'obis', limit = 1, obisopts)


?obis_search()

occ

"basis of record"




dd <- one$obis$data$Lates_niloticus
d <- as.data.frame()
