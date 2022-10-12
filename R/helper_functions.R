# Compute the vector cross product between x and y, and return the components
# indexed by i.
CrossProduct3D <- function(x, y, i=1:3) {
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)

  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1

  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
            x[Index3D(i + 2)] * y[Index3D(i + 1)])
}


#---------------------------------

# Calculate Euclidean distance

#Distance from origin (x1,y1,z1) to insertion (x2,y2,z2); for muscle vectors
#Distance of Muscle (dm) = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)

euc_dist <- function(start, end) {
  sqrt((start[1] - end[1])^2 + (start[2] - end[2])^2 + (start[3] - end[3])^2)
}

#---------------------------------

# Function to calculate crossproduct of 3D force (or unit) vectors
#taken from Matthew Lundberg on StackOverflow (https://stackoverflow.com/questions/36798301/r-compute-cross-product-of-vectors-physics)

xprod <- function(...) {
  args <- list(...)

  # Check for valid arguments

  if (length(args) == 0) {
    stop("No data supplied")
  }
  len <- unique(sapply(args, FUN=length))
  if (length(len) > 1) {
    stop("All vectors must be the same length")
  }
  if (len != length(args) + 1) {
    stop("Must supply N-1 vectors of length N")
  }

  # Compute generalized cross product by taking the determinant of sub-matricies

  m <- do.call(rbind, args)
  sapply(seq(len),
         FUN=function(i) {
           det(m[,-i,drop=FALSE]) * (-1)^(i+1)
         })
}

#------------------------------------

# Calculate magnitude of a 3D vector

vect_mag <- function(x) {
  sqrt(x[1]^2 + x[2]^2 + x[3]^2)
}
