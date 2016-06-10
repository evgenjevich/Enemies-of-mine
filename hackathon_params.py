import json


def hackParams(N, stp, dx):
    
    params = {'N' : N,  'steps' : stp, 'dx' : dx}

    with open('params.json', 'w') as fp:
        json.dump(params, fp)
