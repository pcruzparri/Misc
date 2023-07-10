import os 

os.system("sudo chmod 660 /dev/spi*")
os.system("sudo chown root:spi /dev/spi*")
os.system("sudo chmod 660 /dev/i2c*")
os.system("sudo chown root:i2c /dev/i2c*")
os.system("sudo chmod 660 /dev/gpiomem")
os.system("sudo chown root:gpio /dev/gpiomem")
