from PIL import Image

shape = "bolt"

image = Image.open(f"{shape}.png")
image = image.convert("L")
width, height = image.size

white_pixel_indices = []

for y in range(height):
    for x in range(width):
        pixel_value = image.getpixel((x, y))
        if pixel_value == 255:
            pixel_index = y * width + x
            white_pixel_indices.append(pixel_index)

with open(f"{shape}.pat", "w") as output_file:
    for index in white_pixel_indices:
        output_file.write(f"{index} 1.000000\n")

