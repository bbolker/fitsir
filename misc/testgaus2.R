harbin2 <- setNames(harbin,c("times","count"))
(f1 <- fitsir(harbin2, 
               start=c(beta=2, gamma=1, N=2e3, i0=0.0001, sigma=10),
               type="death"))
plot(f1)

f1_g2 <- fitsir(harbin2, 
              start=c(beta=2, gamma=1, N=2e3, i0=0.0001),
              dist="gaussian2",
               type="death")

f1_g2 <- fitsir(harbin2, 
              start=c(beta=2, gamma=1, N=2e3, i0=0.0001),
               type="death")
