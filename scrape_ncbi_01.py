# from v5 copy
# NCBI scrapping

import requests, time, json, os, traceback
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

ncbi_dataset_file = "ncbi_dataset_final_01.json"

main_link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Tree&id=2153479&lvl=6&keep=1&srchmode=1&unlock'
soup = req_func(main_link)

main_li = soup.find('li', type='circle')
species_li_all = main_li.find_all('li', type='square')

if os.path.isfile(ncbi_dataset_file):
    print('exsisting file found')
    with open(ncbi_dataset_file) as file:
        data = json.load(file)
        last_index = data['species'][-1]['index'] + 1
        ncbi_dataset = data
else: 
    print('no file found, starting from index 0')
    last_index = 0 
    ncbi_dataset = {"species": []}


for i in  range (last_index, len(species_li_all)):
    a_tag = species_li_all[i].find('a')

    species_link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/' + a_tag['href']
    print(species_link, end='\n\n')
    
    try:
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

    except Exception as e:
        print(traceback.format_exc(), end='\n\n')
        with open(ncbi_dataset_file, 'w') as file:
            json.dump(ncbi_dataset, file, indent=2)
        
        print("\n[error at: %s s]" % (time.time() - start_time))
        
        quit()

with open(ncbi_dataset_file, 'w') as file:
            json.dump(ncbi_dataset, file, indent=2)

# insertion_sorting
with open(ncbi_dataset_file) as json_file:
  a = json.load(json_file)
  list_a = a['species']
  print(len(list_a)) # 15798

sorted_ncbi_dataset = {}
sorted_ncbi_dataset.setdefault("species", insertion_sort(list_a))
print(len(insertion_sort(list_a)))

with open('ncbi_dataset_sorted_final_01.json', 'w') as file:
        json.dump(sorted_ncbi_dataset, file, indent=2)

print("\n[Finished in %s s]" % (time.time() - start_time))