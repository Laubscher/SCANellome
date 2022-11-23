import base64
with open("blue.png", "rb") as image2string:
  converted_string = base64.b64encode(image2string.read())
print(converted_string)
