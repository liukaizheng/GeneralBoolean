import json
import os

with open ('./compile_commands.json', 'r') as f:
    load_list = json.load(f)

main_dict = load_list[0]
main_file_name = main_dict['file']
par_dir = os.path.dirname(main_file_name) + '/'
par_len = len(par_dir)
base_name = main_file_name[par_len:]

find_str = '.o -c '
command_find_str = base_name + find_str
main_command = main_dict['command']
command_index = main_command.index(command_find_str)

extrude_file = [d['file'][par_len:] for d in load_list]
head_extension = [ '.h', '.hxx', '.hpp', '.hh' ]

dir_path = [par_dir]
head_files = []
while dir_path:
    cur_dir = dir_path[0]
    dir_path.remove(cur_dir)
    cur_file_and_dir = [os.path.join(cur_dir, fd) for fd in os.listdir(cur_dir)]
    new_dir = [d for d in cur_file_and_dir if os.path.isdir(d)]
    new_file = []
    for fd in cur_file_and_dir:
        if fd not in new_dir:
            rf = fd[par_len:]
            if rf not in extrude_file and os.path.splitext(rf)[1] in head_extension:
                new_file.append(rf)
    head_files += new_file

    if(cur_dir == par_dir):
        new_dir.remove(os.path.join(cur_dir, 'build'))

    dir_path += new_dir


for f in head_files:
    nf = os.path.join(par_dir, f)
    file_command = main_command[:command_index] + f + find_str + nf
    load_list.append({
        'directory': main_dict['directory'],
        'command': file_command,
        'file': nf})

with open('compile_commands.json', 'w') as f:
    json.dump(load_list, f)
