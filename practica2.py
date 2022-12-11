import urllib.request
import ssl
import regex as re

''' FUNCIONES '''


def genera_diccionario_enzimas(link):
    # Saltamos las 10 primeras líneas del archivo
    comienzo = 0

    while comienzo < 10:
        next(link)
        comienzo += 1
    # ER
    expresion_nombre = r'([A-Z]([a-z]|[A-Z]|\d)+)(?=\s)'
    expresion_diana = r'((?<=\s)([A-Z]|\^)+)'
    expresion_r = "R"
    expresion_y = "Y"
    expresion_m = "M"
    expresion_k = "K"
    expresion_s = "S"
    expresion_w = "W"
    expresion_b = "B"
    expresion_d = "D"
    expresion_h = "H"
    expresion_v = "V"
    expresion_n = "N"
    er_simbolo = re.compile(r'\^')
    er_nombre = re.compile(expresion_nombre)
    er_diana = re.compile(expresion_diana)
    er_r = re.compile(expresion_r)
    er_y = re.compile(expresion_y)
    er_m = re.compile(expresion_m)
    er_k = re.compile(expresion_k)
    er_s = re.compile(expresion_s)
    er_w = re.compile(expresion_w)
    er_b = re.compile(expresion_b)
    er_d = re.compile(expresion_d)
    er_h = re.compile(expresion_h)
    er_v = re.compile(expresion_v)
    er_n = re.compile(expresion_n)

    diccionario = {}
    nombre = ""
    diana = ""

    # Leemos a partir de la línea 10
    for linea in link:
        posiciones = []
        posicion = 0
        # Quitamos caracteres de inicio y fin
        linea = linea.decode().strip()

        # Obtenemos el nombre
        r = er_nombre.match(linea)
        if r:
            nombre = r.group(0)
        # Obtenemos la diana
        r = er_diana.search(linea)
        if r:
            diana = r.group(0)

        # Posicion de ^
        m = er_simbolo.search(diana)
        if m :
            posicion = m.start()   
            diana = er_simbolo.sub("", diana)    
               
        # Reemplazamos codigo extendido de ADN
        diana = er_r.sub(r'[AG]', diana)
        diana = er_y.sub(r'[CT]', diana)
        diana = er_m.sub(r'[AC]', diana)
        diana = er_k.sub(r'[GT]', diana)
        diana = er_s.sub(r'[CG]', diana)
        diana = er_w.sub(r'[AT]', diana)
        diana = er_b.sub(r'[CGT]', diana)
        diana = er_d.sub(r'[AGT]', diana)
        diana = er_h.sub(r'[ACT]', diana)
        diana = er_v.sub(r'[ACG]', diana)
        diana = er_n.sub(r'[ACGT]', diana)

        # Comprobamos si la la clave ya esta insertada o no
        if nombre in diccionario:
            viejo = diccionario.pop(nombre)
            if '|' in viejo[0]:
                diana =  viejo[0] + r"|(" + diana + r")"
            else:
                diana = r"(" + viejo[0] + r")|(" + diana + r")"
            posiciones = viejo[1]
        
        posiciones.append(posicion)
        diccionario[nombre] = [diana, posiciones.copy()]
        posiciones.clear()

    # Compilamos las exprexiones regulares y guardamos en el diciconario
    for key in diccionario:
        er = re.compile(list(diccionario[key])[0])
        diccionario[key] = [er, list(diccionario[key])[1]]
    return diccionario


def genera_diccionario_dna(link):
    # Creamos la expresión regular
    enzima = r'(?<=>)[A-Z]\.([A-Z]|[a-z]|\d)+'
    er_enzima = re.compile(enzima)
    espacio = r'\s'
    er_espacio = re.compile(espacio)
    nucleotidos = r'(?<=\s)\d+'
    er_nucleotidos = re.compile(nucleotidos)

    datos = []
    nucleos = ""
    name_gen = ""
    gens = ""
    diccionario = {}

    for linea in link:
        # Quitamos caracteres de inicio y fin
        linea = linea.decode().strip()

        # Obtenemos la enzima
        r = er_enzima.search(linea)
        if r:
            name_gen = r.group(0)

            # Obtenemos el número de nucleótidos
            n = er_nucleotidos.search(linea)
            if n:
                nucleos = n.group(0)
                datos.append(nucleos)
                nucleos = ""
        else:
            if len(linea) == 0:
                # Hemos terminado de leer la lista
                if len(name_gen) > 0:
                    # Guardamos los genes leídos
                    datos.append(gens)
                    diccionario[name_gen] = datos.copy()
                # Reseteamos
                name_gen = ""
                gens = ""
                datos.clear()
            else:
                # Quitamos el espacio y concatenamos
                gens += er_espacio.sub("", linea)
    return diccionario


if __name__ == '__main__':
    '''EJERCICIO 1'''
    print("==============")
    print("Cargando bionet...")
    genes = {}
    enzimas = {}

    try:
        # Descargamos los recursos
        url = 'https://aulavirtual.um.es/access/content/group/1896_G_2022_N_N/PRACTICAS/PRACTICA%202/link_bionet.txt'
        context = ssl.create_default_context()
        context.set_ciphers("DEFAULT")
        link = urllib.request.urlopen(url, context=context)
        enzimas = genera_diccionario_enzimas(link)
    except IOError as e:
        print('link_bionet no disponible:', e)

    ''' EJERCICIO 2 '''
    print("Cargando All_C_genes_DNA.txt...")
    try:
        # Descargamos los recursos
        url = 'https://aulavirtual.um.es/access/content/group/1896_G_2022_N_N/PRACTICAS/PRACTICA%202/All_C_genes_DNA.txt'
        context = ssl.create_default_context()
        context.set_ciphers("DEFAULT")
        link = urllib.request.urlopen(url, context=context)
        genes = genera_diccionario_dna(link)

    except IOError as e:
        print('link_bionet no disponible:', e)
    print("--------------")
    '''
    if 'AbaAI' in enzimas:
        print("está")
        exit(0)
    '''
    '''EJERCICIO 3'''
    while True:
        # 3.a) Pedimos gen por teclado
        nombre_gen = input("Gen >> ")

        # 3.b) Si se introduce la cadena vacía, finalizamos
        if len(nombre_gen) == 0:
            print("==============")
            break
        else:
            if nombre_gen in genes.keys():
                # 3.d) Mostramos la cadena de ADN del gen
                print("-------------- " + list(genes[nombre_gen])[0] + " nucleótidos")
                print(list(genes[nombre_gen])[1])
                print("--------------")

                while True:
                    # 3.d.1) Pedimos una enzima de reconocimiento
                    nombre_enzima = input("Enzima >> ")
                    # 3.d.2) Si se introduce la cadena vacía, finalizamos
                    if len(nombre_enzima) == 0:
                        break
                    else:
                        expresion_regular = False
                        enzima_expresion = r'([A-Z]|[a-z]|[0-9])+'
                        er_enzima = re.compile(enzima_expresion)
                        if not er_enzima.fullmatch(nombre_enzima):
                            # 3.d.4) Tratamos el nombre de la encima como una expresión regular
                            nombre_enzima = re.compile(nombre_enzima)
                            expresion_regular = True
                        encuentra = False
                        for key in enzimas:
                            # 3.d.3) Si la cadena coincide con el nombre de una enzima, imprimimos el mapa
                            # de dianas de dicha enzima en la cadena de ADN anterior
                            if expresion_regular:
                                r = nombre_enzima.fullmatch(key)
                            else:
                                if nombre_enzima == key: r = True
                                else: r = False
                            if r:
                                encuentra = True
                                lista_pos = []
                                for m in enzimas[key][0].finditer(genes[nombre_gen][1]):
                                    if m:
                                        # Si solo hay una diana
                                        if m.lastindex is None:
                                            corte = m.start() + enzimas[key][1][0]
                                        # Si tenemos más de una diana -> hay grupos
                                        else:
                                            corte = m.start() + enzimas[key][1][m.lastindex - 1]
                                            '''
                                            print(m.start())
                                            print(m.lastindex)
                                            print(enzimas[key][1])
'''
                                        lista_pos.append(corte)

                                if len(lista_pos) > 0 or not expresion_regular:
                                    print(key + ' # ', end="")
                                    print(lista_pos)

                        if not encuentra:
                            print("Nombre de enzima incorrecto")

                        print("--------------")
                print("--------------")
            else:
                # 3.c) Informamos al user de que no está el gen
                print("El nombre del gen no se encuentra en el diccionario")
                print("--------------")
