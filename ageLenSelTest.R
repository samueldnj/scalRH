# Test to see if integrated growth selectivity model
# matches up

options( digits = 5 )

Linf  <- 75
K     <- 0.3
CV_L  <- 0.2
t0    <- -0.3
A     <- 20
M     <- 0.2
F     <- 0.1

lenSel50 <- 20
lenSel95 <- 65

# OK, make mean length-at-age curve
ages <- 1:A
lenAge <- Linf * (1 - exp(-K * (ages - t0) ) )

# Now make length bins
lWid  <- 5
lMids <- seq(15,125, by = lWid)


# Make probability mtx
probLenAtAge <- matrix(0, nrow = A, ncol = length(lMids) )

# Now populate
for( a in 1:A )
{
  for( l in 1:length(lMids) )
  {
    lHi <- lMids[l] + lWid/2
    if( l == 1 ) probLenAtAge[a,l] <- pnorm(lHi, lenAge[a], lenAge[a]*CV_L )
    if( l > 1 & l < length(lMids) )
      probLenAtAge[a,l] <- pnorm(lHi, lenAge[a], lenAge[a]*CV_L ) - pnorm(lHi - lWid, lenAge[a], lenAge[a]*CV_L )
    if( l == length(lMids) ) 
      probLenAtAge[a,l] <- 1 - pnorm(lHi - lWid, lenAge[a], lenAge[a]*CV_L)
  }
}

par( mfrow = c(2,2) )

# Now let's plot age/length distribution
plot( x = c(0,A), y = range(lMids), type = "n", axes = F )
  for( a in 1:A )
    points( x = rep(a,length(lMids)), y = lMids, cex = 2*sqrt(probLenAtAge[a,]),
            col = "grey50", pch = 16 )
  lines( x = 1:A, y = lenAge, lty = 2, lwd = 3  )
  axis(side = 1)
  axis(side = 2)

# Now, make the selectivity function
sel_len <- 1 / (1 + exp(-log(19) * (lMids - lenSel50) / (lenSel95 - lenSel50)))

sel_age <- numeric(length = A)
Nat <- numeric(length = A)

for( a in 1:A ) 
{
  sel_age[a] <- sum( probLenAtAge[a,] * sel_len )
  if(a == 1) Nat[a] <- 2000
  if( a > 1 ) Nat[a] <- Nat[a-1] * exp( - M - sel_age[a] * F )
}

# Nat <- rep(1,A)

probHarvAge <- sel_age * Nat / sum( Nat )

probHarvLength <- numeric(length(lMids))
probLenBin <- probHarvLength


Nlt <- numeric( length = length(lMids) )
selectedAge <- numeric(length(lMids))
sel_lenByAge <- numeric(length(lMids))
for( l in 1:length(lMids) )
{
  Nlt[l] <- sum( probLenAtAge[,l] * Nat )
  sel_lenByAge[l] <- sum( probLenAtAge[,l] * sel_age / sum(probLenAtAge[,l]) )
  selectedAge[l] <- sum( probLenAtAge[,l] * Nat *  sel_age )
  probHarvLength[l] <- sum( probLenAtAge[,l] * probHarvAge )
} 

probLenBin <- probHarvLength / sum(probHarvLength)

selectedLength <- sel_len*Nlt

ul <- selectedLength / sum( selectedLength)

plot( selectedLength, lMids, col = 1:length(lMids) )
plot( lMids,selectedAge, col = 1:length(lMids) )

# Now compare the two curves directly
plot( selectedAge,selectedLength, col = 1:length(lMids) )
  abline(a = 0, b = 1)

dev.new()
plot( ul, probLenBin )
  abline(a = 0, b = 1)


