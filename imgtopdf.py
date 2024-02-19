from PIL import Image  # install by > python3 -m pip install --upgrade Pillow  # ref. https://pillow.readthedocs.io/en/latest/installation.html#basic-installation

images = [
    Image.open("/Users/peter/Downloads/" + f)
    for f in ["IMG_6030.jpg"]
]

pdf_path = "/Users/peter/source/Courses/CS532/Activity8/Assignment2handwork.pdf"
    
images[0].save(
    pdf_path, "PDF" ,resolution=100.0, save_all=True, append_images=images[1:]
)
