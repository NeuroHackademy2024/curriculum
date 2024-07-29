
room_list_online = [
    'Alder Auditorium',
    'Alder Commons',
    'Alder 103',
    'Alder 105',
    'Alder 106',
    'Alder 107',
    'Office Hours',
    'Miscellaneous Space']
room_list_inperson = [
    'Alder Auditorium',
    'Alder Commons',
    'Alder 103',
    'Alder 105',
    'Alder 106',
    'Alder 107',
    'Alder Courtyard',
    'Coffe Table']

def get_rooms(n=5, inperson=False):
    import numpy as np
    if inperson:
        room_list = room_list_inperson
    else:
        room_list = room_list_online
    sel = np.random.choice(room_list, n, replace=False)
    for (ii,room) in enumerate(sel):
        print(f"{ii+1}. {room}")
    return None


