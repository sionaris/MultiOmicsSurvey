# Create logo mix using magick
library(magick)
library(png)

# Scale the R logo to 200x200
r_img <- image_read(file.path(home, "Resources/r.png")) |> 
  image_scale("200x200!")

# Scale the Python logo to 200x200
py_img <- image_read(file.path(home, "Resources/python.png")) |> 
  image_scale("200x200!")

# Write them back to disk so we can read them as arrays
image_write(r_img, file.path(home, "Resources/r_200x200.png"))
image_write(py_img, file.path(home, "Resources/python_200x200.png"))

# Read each 200x200 image as a numeric array:  [height, width, channels]
# Typically RGBA => a 4-channel array
r_array   <- readPNG(file.path(home, "Resources/r_200x200.png"))      # shape: 200 x 200 x 4
py_array  <- readPNG(file.path(home, "Resources/python_200x200.png"))  # shape: 200 x 200 x 4

nr <- dim(r_array)[1]  # 200
nc <- dim(r_array)[2]  # 200
# Create a blank result array, same shape
res_array <- array(0, dim = c(nr, nc, 4))

# We want:
#  - the "upper triangle" (row < col) to come from the R logo
#  - the "lower triangle" (row > col) to come from the Python logo
#  - the main diagonal (row == col) to be white

for(i in seq_len(nr)) {
  for(j in seq_len(nc)) {
    
    if(j > i) {
      # Above diagonal => pick from R
      res_array[i, j, ] <- r_array[i, j, ]
      
    } else if(j < i) {
      # Below diagonal => pick from Python
      res_array[i, j, ] <- py_array[i, j, ]
      
    } else {
      # On the diagonal => white background
      # RGBA for white = c(1,1,1,1)
      res_array[i, j, ] <- c(1,1,1,1)
    }
  }
}

# Write out the combined image
writePNG(res_array, file.path(home, "Resources/r_python_mashup.png"))
