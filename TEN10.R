ls();
Lx = 10;
Ly = 10;
N = 3500;

minx = 0;
miny = 0;
maxx = Lx;
maxy = Ly;
mino = -90;
maxo = 90;
mino = mino*pi/180;
maxo = maxo*pi/180;

num <- 1:50000;

for (i in num){
  file <- paste("D:/Data/","TEXT_",toString(i),".txt", sep = "", collapse = NULL)
  Xc = runif(N,minx,maxx);
  Yc <- runif(N,miny,maxy);
  Orientation <- runif(N,mino,maxo);

  Tabl <- data.frame(1:N, Xc, Yc , Orientation);

  write.table(Tabl, file, append = FALSE, sep = "\t\t\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE);
}