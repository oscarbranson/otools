import urllib2
from bs4 import BeautifulSoup
import re
from multiprocessing import Pool
import pandas as pd


# Scrape all the element data from periodictable.com
def elementScrape(path="./"):
    url = "http://www.periodictable.com/Properties/A/StableIsotopes.html"
    urlbase = "http://www.periodictable.com/Properties/A/"

    page = urllib2.urlopen(url)
    soup = BeautifulSoup(page)

    links = soup.find_all("a", href=re.compile('\.\.\/\.\.\/Isotopes'))
    links = [elem['href'] for elem in links]

    # print int(re.search(r"\.\.\/\.\.\/Isotopes\/([0-9]+)\.[0-9]\/index\.html", links[0]).group(1))

    def elscrape(i):
        suburl = urlbase + i
        subpage = urllib2.urlopen(suburl).read()
        subsoup = BeautifulSoup(subpage.decode('utf-8', 'ignore'))

        element = subsoup.br.find_all("table")[3].big.big.big.text
        mass = subsoup.br.findAll('td', align='left', valign="top")[0].text
        abundance = subsoup.br.findAll('td', align='left', valign="top")[1].text[:-1]
        abundance = re.sub(u"1.\\xd7102", "100", abundance)
        electrons = re.search(r"\.\.\/\.\.\/Isotopes\/([0-9]+)\.[0-9]+\/index\.html", i).group(1)
        return (str(element), float(mass), float(abundance), int(electrons))

    pool = Pool(4)
    isotopes = pool.map(elscrape, links)
    pool.close()
    pool.join()

    iso = pd.DataFrame(isotopes, columns=['isotope', 'mass', 'abundance', 'e'])

    iso['element'] = iso['isotope'].str.replace(r'[0-9]', '')
    iso['mass_n'] = iso['isotope'].str.replace(r'[A-Za-z]', '')

    iso.to_pickle(path+'elements.pkl')
