"""
	Created by ilseop, Lee on 2018. 8. 28..
"""

from logging import error, info
from datetime import datetime

def infolog(message):
	info("{0} {1}".format(datetime.now().isoformat(), message))
	print("{0} {1}".format(datetime.now().isoformat(), message))

def errorlog(message):
	error("{0} {1}".format(datetime.now().isoformat(), message))	
	print("{0} {1}".format(datetime.now().isoformat(), message))