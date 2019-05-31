import xmltodict
from dicttoxml import dicttoxml
from xml.dom.minidom import parseString


class XmlHandler:
    @staticmethod
    def str_to_bool(str):
        if str.lower() == 'true':
            return True
        elif str.lower() == 'false':
            return False
        else:
            print('str is wrong')
            return False

    @staticmethod
    def load_param(fpath):
        with open(fpath, 'rt') as f:
            xml = f.read()
        dict = xmltodict.parse(xml)['root']
        if 'bool' in dict:
            for key, value in dict['bool'].items():
                value = XmlHandler.str_to_bool(value)
                dict['bool'][key] = value
        if 'numeric' in dict:
            for key, value in dict['numeric'].items():
                value = float(value)
                dict['numeric'][key] = value
        return dict

    @staticmethod
    def to_xml(dict):
        xml = dicttoxml(dict, attr_type=False)
        with open('parameter.xml', 'wt') as f:
            f.write(parseString(xml).toprettyxml())

    @staticmethod
    def append_day_trade(user_param, symbol, date):
        if user_param['day-trade'] is None:
            user_param['day-trade'] = {symbol: str(date)}
        else:
            user_param['day-trade'][symbol] = str(date)
        return user_param