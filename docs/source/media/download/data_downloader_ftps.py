import ftplib
import getpass

ftps = ftplib.FTP_TLS('archive.stsci.edu')

#https://archive.stsci.edu/ftp.html
#If the data you are downloading is not exclusive access (that is, is publicly accessible), 
#then you can use an anonymous connection, with "anonymous" as the username and your email address as the password. 
#(We don't do anything with the email address, so you can enter something else if you like.)

user = getpass.getpass(prompt='Enter email address (if anonymous -> enter "anonymous"): ', stream=None) 
passwd = getpass.getpass(prompt='Enter password (if anonymous -> enter email address): ', stream=None) 

ftps.login(user=user, passwd=passwd)
ftps.prot_p() # This is a really good idea :)
ftps.cwd('stage')
print(ftps.dir())

ftps.cwd('anonymous/anonymous12345') # can be found in the email

filenames = ftps.nlst()

print(filenames)
for filename in filenames:
    print("getting " + filename)
    with open(filename, 'wb') as fp: 
        ftps.retrbinary('RETR {}'.format(filename), fp.write)

