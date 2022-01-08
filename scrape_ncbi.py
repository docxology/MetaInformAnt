# from v5
# NCBI scrapping

import requests, time, json
from bs4 import BeautifulSoup


def req_func(link):
    r = requests.get(link, timeout=60)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, 'lxml')
    
    return soup

def get_data(small_value= 0):
    fieldset_tag = soup_li.find('fieldset')
    s_name = str(fieldset_tag.find('strong').text.strip())

    small_tag = soup_li.find_all('small')
    if 'txid' in small_tag[small_value].text:
        tax_id_text = small_tag[small_value].text
    else:
        tax_id_text = small_tag[small_value+1].text
    
    tax_id_index = tax_id_text.find('txid')
    tax_id = int(tax_id_text[tax_id_index : ][4 : -1])

    return s_name, tax_id

def insertion_sort(list_a):
    indexing_length = range(1, len(list_a))
    for i in indexing_length:
        value_to_sort = list_a[i]['name']

        while list_a[i-1]['name'] > value_to_sort and i>0:
            list_a[i], list_a[i-1] = list_a[i-1], list_a[i]
            i = i - 1
    return list_a



# main

start_time = time.time()

main_link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Tree&id=2153479&lvl=6&keep=1&srchmode=1&unlock'
soup = req_func(main_link)

main_li = soup.find('li', type='circle')
species_li_all = main_li.find_all('li', type='square')

ncbi_dataset = {"species": []}

for i, x in enumerate(species_li_all):
    a_tag = x.find('a')
    j = i
    
    species_link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/' + a_tag['href']
    print(species_link, end='\n\n')
    
    soup_li = req_func(species_link)
    
    s_name, tax_id = get_data()

    print(f'index serial: {i}')
    print(f'name: {s_name}')
    print(f'txid: {tax_id}')
    print(f'NCBI_link: {species_link}', end='\n\n')

    species_entry ={}

    species_entry.setdefault("name", s_name)
    species_entry.setdefault("index", i)
    species_entry.setdefault("txid", tax_id)
    species_entry.setdefault("ncbi_link", species_link)

    ncbi_dataset['species'].append(species_entry)



# insertion_sorting
sorted_ncbi_dataset = {}

list_a = ncbi_dataset['species']
sorted_ncbi_dataset.setdefault("species", insertion_sort(list_a))

file_name = 'ncbi_dataset_sorted_final.json'
with open(file_name, 'w') as file:
        json.dump(sorted_ncbi_dataset, file, indent=2)


print('total species scraped: ' + str(len(sorted_ncbi_dataset['species'])))
print(f'{file_name} is created on current working directory')
print("\n[Finished in %s s]" % (time.time() - start_time))