import requests
from bs4 import BeautifulSoup


url = 'http://bioengineering.stanford.edu/faculty/'
source_code = requests.get(url)
plain_text = source_code.text
soup = BeautifulSoup(plain_text)
for link in soup.findAll('a', ):
   href = link.get('href')
   print(href)


        


