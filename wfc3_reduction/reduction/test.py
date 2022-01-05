import os
print('/'.join(os.path.realpath(__file__).split('/')))
print('/'.join(os.path.realpath(__file__).split('/')[:-1]))
print('/'.join(os.path.realpath(__file__).split('/')[:-2]))
print('/'.join(os.path.realpath(__file__).split('/')[:-3]))
