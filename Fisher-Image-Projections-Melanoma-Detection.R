library("OpenImageR")
library("plotly")
library("reshape2")
library("pracma")
library("OpenImageR")
library("stats")
library("doParallel")
library("foreach")

setwd("C:/Programming in R")

#--------
#2D Plots
#--------

myImage = "melanoma1.jpg"
original_image = readImage(myImage)
gray_image = rgb_2gray(original_image)

#ZCA_whitening() to whiten the image
?ZCAwhiten()
white_image = ZCAwhiten(gray_image, k = nrow(gray_image), epsilon = 1e-5)

#Display the image
?imageShow
imageShow(myImage)
imageShow(gray_image)
imageShow(white_image)

#Fisher Projection

fisherProjection = function(theta, imageMatrix)
{
  #Rotate and sum columns
  rotatedMatrix = rotateImage(imageMatrix, angle = theta, method = "bilinear")
  projectedData = colSums(rotatedMatrix, na.rm = TRUE)
  
  #Calculate k means
  kiloMeans = kmeans(projectedData, centers = 2, nstart = 50)
  
  #Seperate the projected data into two classes
  pixelData1 = projectedData[kiloMeans$cluster == 1]
  pixelData2 = projectedData[kiloMeans$cluster == 2]
  
  #Return zero if the variance is too small
  if (length(pixelData1) < 2 || length(pixelData2) < 2) 
  {
    return(0)
  }
  
  #Calculate means and variances
  mean1 = mean(pixelData1)
  mean2 = mean(pixelData2)
  var1 = var(pixelData1)
  var2 = var(pixelData2)
  
  #Fisher linear discriminant (μ1 - μ2)^2/(var1 + var2)
  fisherIndex = (mean1 - mean2)^2/(var1 + var2)
}

#Optimal Projection Angle
optimalResult = optimize(
  f = function(theta) fisherProjection(theta, white_image), interval = c(0, 180), maximum = TRUE)

optimalAngle = optimalResult$maximum
paste("Optimal Angle:", optimalAngle)

#Using Plotly to Visualize
?sapply #sapply(list, function)
angles = seq(0, 180, by =2)
fisherValues = sapply(angles, function(theta) fisherProjection(theta, white_image))
paste(fisherValues)

figure = plot_ly(x = ~angles, y = ~fisherValues, type = 'scatter', mode = 'lines', name = 'Fisher Index') %>%
  add_trace(x = ~optimalAngle, y = ~optimalResult$objective, type = 'scatter', mode = 'markers',
            marker = list(color = 'green', size = 15), name = 'Maximum') %>%
  layout(title = 'Optimization Surface for Fisher Index',
         xaxis = list(title = 'Projection Angle (degrees)'),
         yaxis = list(title = 'Fisher Index'))

print(figure)

#Projection and Histogram
rotatedOptimalImage = rotateImage(white_image, optimalAngle, method = "bilinear")
print(rotatedOptimalImage)

optimalProjectionData = colSums(rotatedOptimalImage, na.rm = TRUE)
print(optimalProjectionData)

finalClusters = kmeans(optimalProjectionData, centers = 2)
print(finalClusters)

threshold = mean(finalClusters$centers)
print(threshold)

imageShow(rotatedOptimalImage) #Optimal Image
hist(optimalProjectionData, breaks = 40, col = "red",main = "Histogram of Optimal Projection", xlab = "Projected Values")

#-------
#3D Plot
#-------

smallImage = resizeImage(original_image, width = 64, height = 64)

rgbMatrix = matrix(c(as.vector(smallImage[,,1]),
                     as.vector(smallImage[,,2]),
                     as.vector(smallImage[,,3])),
                   ncol = 3)

almostWhiteMatrix = ZCAwhiten(rgbMatrix, k = 3, epsilon = 1e-5)
whiteMatrix = scale(almostWhiteMatrix, center = FALSE, scale = TRUE)

covMatrix = cov(whiteMatrix)
print(round(covMatrix, 4))


#Fisher Projection
fisherProjection3D = function(phi, theta, whiteData){
  phiRad = phi * pi/180;
  thetaRad = theta * pi/180
  directionVector = c(sin(thetaRad) * cos(phiRad),
                      sin(thetaRad) * sin(phiRad),
                      cos(thetaRad))
  
  projectedData = as.matrix(whiteData) %*% directionVector
  clusters = kmeans(projectedData, centers = 2, nstart = 5)
  class1 = projectedData[clusters$cluster == 1]
  class2 = projectedData[clusters$cluster == 2]
  
  if (length(class1) < 2 || length(class2) < 2) return (0)
  fisherIndex = (mean(class1) - mean(class2))^2 / (var(class1) + var(class2))
  return(fisherIndex)
}

#Optimal 3D Angle
phiAngles = seq(0, 180, by = 10)
thetaAngles = seq(0, 360, by = 10)
fisherMatrix = matrix(NA, nrow = length(phiAngles), ncol = length(thetaAngles))

for (i in 1:length(phiAngles)) {
  for (j in 1:length(thetaAngles)) {
    fisherMatrix[i, j] = fisherProjection3D(phiAngles[i], thetaAngles[j], whiteMatrix)
  }
}

#Optimize
maxFisher = max(fisherMatrix, na.rm = TRUE)
maxIndex = which(fisherMatrix == maxFisher, arr.ind = TRUE)
optimalPhi = phiAngles[maxIndex[1, 1]]
optimalTheta = thetaAngles[maxIndex[1, 2]]

#Viusalize with Plotly (3D Linear Subspace)

figure = layout(
  plot_ly(x = ~thetaAngles, y = ~phiAngles, z = ~fisherMatrix, type = "surface"),
  title = "Fisher Index",
  scene = list(
    xaxis = list(title = "Theta"),
    yaxis = list(title = "Phi"),
    zaxis = list(title = "Fisher Index")
  )
)  

print(figure) #3d Linear Subspace of Fisher Index

#Viusalize with Plotly (Projected Image)

phiRad = optimalPhi * pi/180;
thetaRad = optimalTheta * pi/180
optimalDirectionVector = c(sin(thetaRad) * cos(phiRad),
                           sin(thetaRad) * sin(phiRad),
                           cos(thetaRad)) #Re-use the function one time :D

optimalProjectedData = as.matrix(whiteMatrix) %*% optimalDirectionVector
projectedImage = matrix(optimalProjectedData, nrow = nrow(smallImage), ncol = ncol(smallImage))

finalClusters = kmeans(optimalProjectedData, centers = 2)
threshold = mean(finalClusters$centers)

imageShow(projectedImage) #Projected image

#Viusalize with Plotly (Histogram)
hist(optimalProjectedData, breaks = 40, col = "pink",main = "Histogram of Optimal Projection", xlab = "Projected Values")

#Orthogonal Vectors (Part 1: Copy down the optimal vector)
phiRad = optimalPhi * pi/180;
thetaRad = optimalTheta * pi/180
vector1 = c(sin(thetaRad) * cos(phiRad),
            sin(thetaRad) * sin(phiRad),
            cos(thetaRad)) #Already calculated earlier

#Orthogonal Vector (Part 2: Change of Basis aka: orthogonal to the old vector)

potentialVector = if (abs(vector1[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0) #We dont want the vector to straight up and down
vectorA = cross(vector1, potentialVector)
vectorA = vectorA/sqrt(sum(vectorA^2))

vectorB = cross(vector1, vectorA)
vectorB = vectorB/sqrt(sum(vectorB^2))

#Fisher function
fisherProjectionOrthogonal <- function(alpha, whiteData, basisA, basisB) {
  alphaRad = alpha * pi / 180
  
  #Look for another orthognal vector
  vector2 = cos(alphaRad)*basisA + sin(alphaRad)*basisB
  
  projectedData = as.matrix(whiteData) %*% vector2
  clusters = kmeans(projectedData, centers = 2, nstart = 5)
  class1 = projectedData[clusters$cluster == 1]
  class2 = projectedData[clusters$cluster == 2]
  
  if (length(class1) < 2 || length(class2) < 2) return (0)
  fisherIndex = (mean(class1) - mean(class2))^2 / (var(class1) + var(class2))
  return(fisherIndex)
}

#Find the orthogonal angle
ortho = optimize(
  f = function(alpha) fisherProjectionOrthogonal(alpha, whiteMatrix, vectorA, vectorB),
  interval = c(0, 360),
  maximum = TRUE
)
optimalAlpha = ortho$maximum
print(optimalAlpha)

#(Orthogonal Vectors part 3) Visualize with plotly the second optimal projection

optimalAlphaRad = optimalAlpha * pi/180
thetaRad = optimalTheta * pi/180
optimalVector2 = cos(optimalAlphaRad)*vectorA + sin(optimalAlphaRad)*vectorB

optimalProjectedDataOrtho = as.matrix(whiteMatrix) %*% optimalVector2
projectedImageOrtho = matrix(optimalProjectedDataOrtho, nrow = nrow(smallImage))

finalClustersOrtho = kmeans(optimalProjectedDataOrtho, centers = 2)
thresholdOrtho = mean(finalClustersOrtho$centers)

imageShow(projectedImageOrtho) #Second Projection
hist(optimalProjectedDataOrtho, breaks = 40, col = "yellow", main = "Histogram 
     of 2nd Projection", xlab = "Projected Values") #Second Histogram

#3D orthogonal plot with plotly
alphaPlotAngles = seq(0, 360, by = 2) #sequence of angles

fisherValuesOrtho =
  sapply(alphaPlotAngles, function(a) 
    fisherProjectionOrthogonal(
      a, whiteMatrix, vectorA, vectorB)) #Fisher for the angles

alphaPlotAnglesRad = alphaPlotAngles*(pi/180)

figureOrtho3D = layout(
  plot_ly(
    x = ~cos(alphaPlotAnglesRad), y = ~sin(alphaPlotAnglesRad), z = ~fisherValuesOrtho,
    type = "scatter3d", mode = "lines+markers", marker = list(size = 4), line = list(width = 2)
  ),
  title = "Fisher Index",
  scene = list(
    xaxis = list(title = "cos(α)"),
    yaxis = list(title = "sin(α)"),
    zaxis = list(title = "Fisher Index")
  )
)

print(figureOrtho3D)

#Third orthogonal projection

optimalVector3 = cross(vector1, optimalVector2)
optimalVector3 = optimalVector3/sqrt(sum(optimalVector3^2))
print(optimalVector3)

optimalProjectedData3 = as.matrix(whiteMatrix) %*% optimalVector3
projectedImage3 = matrix(optimalProjectedData3, nrow = nrow(smallImage),)

finalClusters3 = kmeans(optimalProjectedData3, centers = 2)
threshold3 = mean(finalClusters3$centers)

imageShow(projectedImage3) #Third Projection
hist(optimalProjectedData3, breaks = 40, col = "violet", main = "Histogram 
     of 3rd Projection", xlab = "Projected Values") #Third Histogram



myImage = "melanoma1.jpg"
original_image = readImage(myImage)
original_image = resizeImage(original_image, width = 64, height = 64)


optimalProjection = function(original_image, nstart = 5, niter = 25){
  imageMatrix = matrix(c(as.vector(original_image[,,1]),
                         as.vector(original_image[,,2]),
                         as.vector(original_image[,,3])),
                       ncol = 3)
  almostWhiteMatrix = ZCAwhiten(imageMatrix, k = 3, epsilon = 1e-5)
  whiteMatrix = scale(almostWhiteMatrix, center = FALSE, scale = TRUE)
  
  #Core usage
  cores = detectCores()
  cluster = makeCluster(cores - 1)
  registerDoParallel(cluster)
  
  phiAngles = seq(0, 180, by = 10)
  thetaAngles = seq(0, 360, by = 10)
  angleGrid = expand.grid(phi = phiAngles, theta = thetaAngles)
  
  #angleGrid Search
  fisherValues = foreach(i = 1:nrow(angleGrid), .combine = "c", .export = c("whiteMatrix", "nstart", "niter", "angleGrid")) %dopar%{
    phi = angleGrid$phi[i]
    theta = angleGrid$theta[i]
    
    phiRad = phi*pi/180
    thetaRad = theta*pi/180
    newVector = c(sin(thetaRad)*cos(phiRad),
                  sin(thetaRad)*cos(thetaRad),
                  cos(thetaRad))
    projectionData = as.matrix(whiteMatrix) %*% newVector
    
    #New cluster with Niter
    clusters = kmeans(projectionData, centers = 2, nstart = nstart, iter.max = niter)
    class1 = projectionData[clusters$cluster == 1]
    class2 = projectionData[clusters$cluster == 2]
    
    #Fisher Index
    fisherIndex = (mean(class1) - mean(class2))^2 / (var(class1) + var(class2))
    return(fisherIndex)
  }
  stopCluster(cluster)
  
  maxFisher = max(fisherValues, na.rm = TRUE)
  bestFisher = which.max(fisherValues)
  optimalPhi = angleGrid$phi[bestFisher]
  optimalTheta = angleGrid$theta[bestFisher]
  
  #Optimal Vector 1
  phiRad = optimalPhi * pi/180
  thetaRad = optimalTheta * pi/180
  vector1 = c(sin(thetaRad) * cos(phiRad), sin(thetaRad) * sin(phiRad), cos(thetaRad))
  
  #Optimal vector 2
  newVector = if (abs(vector1[3]) < .9) c(0, 0, 1) else c(1, 0, 0)
  #Normalize Vectors
  vectorA = pracma::cross(vector1, newVector)
  vectorA = vectorA / sqrt(sum(vectorA^2))
  
  vectorB = pracma::cross(vector1, vectorA)
  vectorB = vectorB / sqrt(sum(vectorB^2))
  
  #Orthonal Fisher Function
  
  fisherProjectionOrthogonal = function(alpha, whiteMatrix2, nstart2, niter2, vectorA, vectorB){
    alphaRad = alpha*pi / 180
    vector2 = cos(alphaRad)*vectorA+sin(alphaRad)*vectorB
    projectionData = as.matrix(whiteMatrix2) %*% vector2
    clusters = kmeans(projectionData, centers = 2, nstart = nstart2, iter.max = niter2)
    class1 = projectionData[clusters$cluster == 1]
    class2 = projectionData[clusters$cluster == 2]
    
    fisherIndex = (mean(class1) - mean(class2))^2 / (var(class1) + var(class2))
    return(fisherIndex)
  }
  
  orthogonalResult = optimize(
    f = function(a) fisherProjectionOrthogonal(a, whiteMatrix, nstart, niter, vectorA, vectorB),
    interval = c(0, 360), maximum = TRUE
  )
  optimalAlpha = orthogonalResult$maximum
  optimalAlphaRad = optimalAlpha*pi/180
  
  vector2 = cos(optimalAlphaRad) * vectorA + sin(optimalAlphaRad) * vectorB
  
  #Optimal vector 3
  vector3 = pracma::cross(vector1, vector2)
  vector = vector3/sqrt(sum(vector3^2))
  
  #Generating Images
  projectionData1 = as.matrix(whiteMatrix) %*% vector1
  projectionData2 = as.matrix(whiteMatrix) %*% vector2
  projectionData3 = as.matrix(whiteMatrix) %*% vector3
  
  imgDim = dim(original_image)
  
  image1 = matrix(projectionData1, nrow = imgDim[1], ncol = imgDim[2])
  image2 = matrix(projectionData2, nrow = imgDim[1], ncol = imgDim[2])
  image3 = matrix(projectionData3, nrow = imgDim[1], ncol = imgDim[2])
  
  return(list(proj1 = image1, proj2 = image2, proj3 = image3))
} 

#presenting the images
projectedImages = optimalProjection(original_image)

plot.new(); imageShow(projectedImages$proj1) #Projection 1
plot.new(); imageShow(projectedImages$proj2) #Projection 2
plot.new(); imageShow(projectedImages$proj3) #Projection 3