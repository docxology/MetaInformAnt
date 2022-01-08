# from v7


import requests, time, json
from bs4 import BeautifulSoup


def req_func(link):
    r = requests.get(link, timeout=60)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, 'lxml')
    
    return soup

def rename_keys(target_element):
    keys = list(target_element.keys())
    for key in keys:
        if not key == 'Page':
            Temp_value = target_element[key]
            del target_element[key]
            new_key = table_names[x-1] + '_' + key
            target_element[new_key] = Temp_value

    return target_element

def get_json(x, tr_list):
    tr_link = tr_list[x].find('td', class_='cargo-tablelist-tablename').a['href']
    
    response = requests.get(tr_link, timeout=60)
    soup = BeautifulSoup(response.text, 'lxml')
    
    more_div = soup.find('div', id='mw-content-text')
    
    link = 'https://antwiki.org' + more_div.find('a', title='Special:CargoQuery')['href'] + '&format=json'
    # print(link)
    json_soup = req_func(link)
    view_json_div = json_soup.find('div', id='mw-content-text')
    links = view_json_div.find_all('a')

    for y in links:
        if y.text == 'View JSON':
            json_link = y['href']

    json_link = json_link.replace('limit=100', 'limit=5000')

    response = requests.get(json_link, timeout=60)
    data_pull = response.json()

    return data_pull

def special_datapull(x, tr_list, next_link=False, items=[]):
    if not next_link:
        tr_link = tr_list[x].find('td', class_='cargo-tablelist-tablename').a['href']
            
        response = requests.get(tr_link, timeout=60)
        soup = BeautifulSoup(response.text, 'lxml')

        more_div = soup.find('div', id='mw-content-text')

        link = 'https://antwiki.org' + more_div.find('a', title='Special:CargoQuery')['href']
        link = link.replace('limit=100', 'limit=5000')
        link = link.replace('offset=100', 'offset=0')

    else:
        link = next_link

    response = requests.get(link, timeout=60)
    soup = BeautifulSoup(response.text, 'lxml')

    main_div = soup.find('div', id='mw-content-text')
    main_div = main_div.find('div')
    table_div = soup.find('table', class_='cargoTable noMerge sortable') 
    all_tr = table_div.find('tbody')
    all_tr = all_tr.find_all('tr') #  5000 chunk per loop 

    for tr in all_tr:
        all_td = tr.find_all('td')
        item = {}
        for td in all_td:
            key = str(td['class'])[2:-2].replace('field_', '')
            value = td.text

            item.setdefault(key, value)

        items.append(item)

    next_link = False
    all_p = main_div.find_all('p') # next 5000 link <p>
    for p in all_p:
        if p.find('a', title='Next 5,000 results'):
            next_link = 'https://antwiki.org' + p.find('a', title='Next 5,000 results')['href']

    if next_link:
        special_datapull(x, tr_list, next_link, items)

    return items

def listing_element(element, main_list):
    listing_element = list()
    items = list(element.items())
    for item in items:
        listing_element.append(list(item))
    
    fields_len = len(listing_element)
    if len(main_list) == 0:
         main_list.append(element)
    else:
        trigger = False
        for x in main_list:
            if x[listing_element[0][0]] == listing_element[0][1]:
                trigger = True
                for field in range(1, fields_len):
                    if listing_element[field][0] in list(x.keys()):
                        if type(x[listing_element[field][0]]) != list:
                            temp_val = list()
                            temp_val.append(x[listing_element[field][0]])
                        else:
                            temp_val = x[listing_element[field][0]]
                        x[listing_element[field][0]] = list()
                        for temp in temp_val:
                            x[listing_element[field][0]].append(temp)
                        x[listing_element[field][0]].append(listing_element[field][1])
                    else:
                        x[listing_element[field][0]] = listing_element[field][1]
        if not trigger:
            main_list.append(element)
        
    return main_list

def insertion_sort(list_a):
    indexing_length = range(1, len(list_a))
    for i in indexing_length:
        value_to_sort = list_a[i]['Page']

        while list_a[i-1]['Page'] > value_to_sort and i>0:
            list_a[i], list_a[i-1] = list_a[i-1], list_a[i]
            i = i - 1
    return list_a

# main

start_time = time.time()

antWiki_link = 'https://www.antwiki.org/wiki/Special:CargoTables'
soup = req_func(antWiki_link)

main_table = soup.find('table', class_='mw-datatable cargo-tablelist')
tr_list = main_table.find_all('tr')

antwiki_dataset = {"antwiki": []}

# working tables
table_names = ['Associate', 'Ataglance', 'Economolab3D', 'FlightMonth', 'FossilFormation', 'FossilOccurrence', 'Karyotype', 'MaleMorphology', 'TaxonName', 'TypeSpecimen', 'WorkerMorphology']

for x in range(1, len(tr_list)): # range must start from 1
    
    if x != 5:        
        print('scraping: ' + table_names[x-1])
        if x == 9:
            items_main_09 = []
            data_pull_09 = special_datapull(9, tr_list, False, items_main_09)

        elif x == 10:
            items_main_10 = []
            data_pull_10 = special_datapull(10, tr_list, False, items_main_10)

        else:
            data_pull = get_json(x, tr_list)


        if x != 9 and x != 10:
            for target_element in data_pull:
                target_element = rename_keys(target_element)

                main_list = listing_element(target_element, antwiki_dataset['antwiki'])
                antwiki_dataset['antwiki'] = main_list
        
        elif x == 9:
            for target_element in data_pull_09:
                target_element = rename_keys(target_element)

                main_list = listing_element(target_element, antwiki_dataset['antwiki'])
                antwiki_dataset['antwiki'] = main_list
        
        elif x == 10:
            for target_element in data_pull_10:
                target_element = rename_keys(target_element)

                main_list = listing_element(target_element, antwiki_dataset['antwiki'])
                antwiki_dataset['antwiki'] = main_list

print('antwiki dataset created...')


# insertion sort & adding 'txid': 0 by default
list_a = antwiki_dataset['antwiki']

for page in list_a:
    page.setdefault("txid", 0)

print('insertion sorting starts\n')
sorted_antwiki_dataset = {}
sorted_antwiki_dataset.setdefault("antwiki", insertion_sort(list_a))

file_name = "antwiki_dataset_sorted_final_01.json"
with open(file_name, 'w') as file:
        json.dump(sorted_antwiki_dataset, file, indent=2)

print('sorted antwiki_dataset len --> ' + str(len(sorted_antwiki_dataset['antwiki'])))
print(f'{file_name} is created on current working directory')



print("\n\n[Finished in %s s]" % (time.time() - start_time))