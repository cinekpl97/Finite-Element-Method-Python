import json


def get_settings():
    with open('D:\\PROGRAMOWANIE\\Python programy\\Finite-Element-Method\\structures\\settings.json') as file:
        settings = json.load(file)
    return settings
