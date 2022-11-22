import python_forex_quotes
import json
from datetime import datetime
import MySQLdb


def onUpdate(value):
    try:
        dataset = json.loads(value)
        symbol = dataset['s']
        dataset = {"p": dataset['p'], "a": dataset['a'], "b": dataset['b'], "t": datetime.fromtimestamp(dataset['t']/1000)}
        data_forex[symbol].append(dataset)
    except Exception as exc:
        print(f"Error {exc}")



def onMessage(value):
    print(value)


def onConnect():
    client.subscribeToAll()


if __name__ == "__main__":
    dbconnect = MySQLdb.connect("localhost", "test", "testpwd123", "forex")
    try:
        cursor = dbconnect.cursor()
        cursor.execute("SELECT VERSION()")

        data = cursor.fetchone()
        if data:
            print('Version retrieved: ', data)
        else:
            print('Version not retrieved.')


        client = python_forex_quotes.ForexDataClient("woV3dF1LV9WXFb4gElnoCd76TtP0Oe7e")
        symbols = client.getSymbols()

        script = f"INSERT INTO `instruments` (`symbol`) VALUES"
        for item in symbols:
            script += f"('{item}'),"
        script = script[:-1]
        script += ";"

        try:
            cursor.execute(script)
            dbconnect.commit()
        except Exception as exc:
            print(f"{exc}")
            dbconnect.rollback()

        data_forex = {}
        for item in symbols:
            data_forex[item] = []

        client.onUpdate(onUpdate)
        client.onConnect(onConnect)
        client.onMessage(onMessage)
        client.connect()
    except:
        pass
    finally:
        dbconnect.close()
