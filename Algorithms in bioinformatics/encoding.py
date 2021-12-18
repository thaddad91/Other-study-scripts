#!/usr/bin/env python
"""
This is a script to encode numeric messages.
The messages of five numbers encode a day and a place.
E.g. the message 34681 codes for saturday, fountain.
"""

from random import choice

def encode(day, place):
    """
    Return encrypted message as a list of 5 integers
    day: string, 'monday' through 'sunday'
    place: string, 'fountain', 'tree' or 'church'
    """
    valid = False
    while not valid:
        encoded_day = encode_day(day)
        print("day: ",encoded_day)
        if sum(encoded_day) %2 == 0:
            must_be_even = True
        else:
            must_be_even = False
        encoded_place = encode_place(place, must_be_even)
        print("place: ",encoded_place)
        encoded_message = encoded_day + encoded_place
        valid = is_valid(encoded_message)
    return encoded_message

def is_valid(code):
    """
    Return True if code is valid, False otherwise
    code: tuple or list of five integers i where 0 <= 1 < 10
    This function performs two initial checks:
    If the fourth number equals the fifth number, the code is invalid
    If the sum of all numbers is odd, the code is invalid
    Otherwise, it may be valid
    """
    if code[3] == code[4]:
        return False
    if sum(code)%2 == 0:
        return True
    return False

def encode_day(day):
    """
    Return first 3 numbers of code as list of integers
    day: string, 'monday' through 'sunday'
    Numbers corresponding to days:
    monday: 1, tuesday: 2, ...., sunday: 7
    """
    days = {
        'monday': 1,
        'tuesday': 2,
        'wednesday': 3,
        'thursday': 4,
        'friday': 5,
        'saturday': 6,
        'sunday': 7,
        }
    dnum = days[day]
    first, second = choice(product_options(dnum, dnum+9))
    third = first*second-dnum
    return [first, second, third]

def product_options(min_product, max_product):
    """
    Return options (a,b) such that min_product <= a*b <= max_product
    min_product: integer, lower bound on the product
    max_product: integer, upper bound on the product
    """
    options = []
    for x in range(1,10):
        for y in range(1,10):
            if x*y > max_product:
                break
            if x*y < min_product:
                continue
            options.append((x,y))
    return options

def encode_place(place, must_be_even):
    """
    Return last 2 numbers of code as list of integers
    place: string, 'fountain', 'tree' or 'church'
    must be even. boolean indicating whether the sum of both
    returned numbers should be even.
    numbers corresponding to places:
    fountain: 2, tree: 7, church: 8
    """
    valid_places = {
        'fountain': 2,
        'tree': 7,
        'church': 8,
        }
    place_code = valid_places[place]
    is_ok = False
    while not is_ok:
        chosen = choice([1, -1])
        fifth = place_code + chosen
        if chosen == 1:
            fourth = choice(range(1, fifth))
        else:
            fourth = choice(range(fifth+1, 10))
        if (fourth+fifth) % 2 == 0 and must_be_even:
            is_ok = True
        if (fourth+fifth) % 2 != 0 and must_be_even is False:
            is_ok = True
    return [fourth, fifth]

def main():
    my_code = encode('sunday', 'fountain')
    print("encoded day and place: ", my_code)

main()
