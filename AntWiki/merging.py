import json, time


def binary_search(main_list, item):
    begin_index = 0
    end_index = len(main_list) - 1

    while begin_index <= end_index:
        mid_point = begin_index + (end_index - begin_index) // 2
        mid_point_value = main_list[mid_point]['Page']

        if mid_point_value == item['name']:
            return mid_point
        
        elif item['name'] < mid_point_value:
            end_index = mid_point - 1
        
        elif item['name'] > mid_point_value:
            begin_index = mid_point + 1
        
    return None

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

ncbi_file = "ncbi_dataset_sorted_final_01.json"
with open(ncbi_file) as json_file:
  a = json.load(json_file)
  ncbi_list = a['species']

antwiki_file = "antwiki_dataset_sorted_final_01.json"
with open(antwiki_file) as json_file:
  b = json.load(json_file)
  antwiki_list = b['antwiki']

print(len(antwiki_list))
print(len(ncbi_list), end='\n\n')

not_matched = 0
matched = 0

for i, ncbi_item in enumerate(ncbi_list):
    print('len at begin: ' + str(len(antwiki_list)))
    return_value = binary_search(antwiki_list, ncbi_item)
    
    if return_value == None:
      print('not matched')
      
      temp_value = ncbi_item['name'] 
      del ncbi_item['name']
      ncbi_item['Page'] = temp_value
      antwiki_list.append(ncbi_item)
      not_matched = not_matched + 1
    
    else:
      print('matched')
      antwiki_list[return_value]['txid'] = ncbi_item['txid']
      print(antwiki_list[return_value])
      print(ncbi_item)
      matched = matched + 1

    print(f'{i} cycle ends', end='\n\n')

print('insertion sorting starts\n')
sorted_antwiki_dataset = {}
sorted_antwiki_dataset.setdefault("antwiki", insertion_sort(antwiki_list))

file_name = 'merged_sorted_final_01.json'
with open(file_name, 'w') as file:
        json.dump(sorted_antwiki_dataset, file, indent=2)


print('report:\n')
print('merged antwiki datatest len--> ' + str(len(sorted_antwiki_dataset['antwiki'])))

print(f'\nnot matched in ncbi dataset--> {not_matched}')
print(f'matched in ncbi dataset--> {matched}')
print(f'total: {not_matched + matched}')
print(f'species in ncbi dataset --> {len(ncbi_list)}\n')

print(f'{file_name} is created on current working directory')

print("\n\n[Finished in %s s]" % (time.time() - start_time))