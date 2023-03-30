"""
21849    v2841(-4, 1)
27964    v3284(-4, 3)
30829     v3478(1, 5)
"""

import pandas as pd
import json, os, hashlib, re

dir = '/pkgs/tmp/proofs/'

def save_as_file(row):
    file = open(dir + row['name'], 'w')
    file.write(row['proof_nonord'] + '\n')
    file.close()
    
def save_redis_result(result):
    result = json.loads(result)
    try:
        proof = result['examples'][2][0]
        name = json.loads(proof)['name'].replace(',', ', ')
        file = open(dir + name, 'w')
        file.write(proof + '\n')
        file.close()
        print(name)
    except KeyError:
        pass

def done():
    return os.listdir(dir)

def common_proofs():
    ans = []
    for file in done():
        proof = json.loads(open(dir + file).read())['proof']
        ans.append(hashlib.md5(json.dumps(proof)).hexdigest())
    return ans

def load_proof(file):
    return json.loads(open(dir + file).read())['proof']


def check_labels(file):
    proof = load_proof(file)
    for labels, word in proof:
        labels = labels.split('.')
        if len(set(labels)) != len(labels):
            assert False
    
                
def fix_proof(proof):
    def invert_word(word):
        return word.swapcase()[::-1]

    if isinstance(proof, str):
        proof = json.loads(proof)
    claims = proof['proof']
    claims = [(a.split('.'), b) for a, b in claims]
    fixed = []
    for path, word in claims:
        new_path = []
        for g in path:
            ginv = invert_word(g)
            if ginv in new_path:
                i = new_path.index(ginv)
                new_path = new_path[:i]
            new_path.append(g)
        fixed.append(('.'.join(new_path), word))

    proof['proof'] = fixed
    return json.dumps(proof)

def fix_all_proofs():
    for file in os.listdir(dir):
        original = json.loads(open(dir + file).read())
        corrected = fix_proof(original).replace(' ', '')
        final = re.sub('"name":".*?,', lambda x:x.group().replace(',', ', '), corrected)
        assert final.find(file) > 0
        newfile = open('/pkgs/tmp/proofs_new/' + file, 'w')
        newfile.write(final + '\n')
        newfile.close()
