'''Defines fuctions and data for tests.'''

__all__ = [
        'oscillator_1st_deriv',
        'oscillator_2nd_deriv',
        'oscillator_euler_t',
        'oscillator_euler_x1',
        'oscillator_euler_x2',
          ]


def oscillator_1st_deriv(t, X):
    ''''$\begin{bmatrix}X\end{bmatrix} =
        \begin{bmatrix}x\\v\end{bmatrix}$
        $\dot{\begin{bmatrix}X\end{bmatrix}} =
        \begin{bmatrix}v\\-\omega^2 x\end{bmatrix}$'''
    w = 1.5
    x, v = X
    dx, dv = [
        v,
        -(w**2) * x
            ]
    dX = [dx, dv]
    return dX


def oscillator_2nd_deriv(t, X):
    '''$\begin{bmatrix}X\end{bmatrix} =
        \begin{bmatrix}x\end{bmatrix}$
        $\ddot{\begin{bmatrix}X\end{bmatrix}} =
        \begin{bmatrix}-\omega^2 x\end{bmatrix}$'''
    w = 1.5
    x = X
    ddx = -(w**2) * x
    ddX = [ddx]
    return ddX


oscillator_euler_t = [
    0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
    1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1,
    3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
    4.7, 4.8, 4.9, 5]
oscillator_euler_x1 = [
    0, 0.1, 0.2, 0.29775, 0.391, 0.477550625, 0.55530375, 0.6223119859375,
    0.6768258875, 0.717337769378906, 0.742621068789062, 0.751764268388193,
    0.74419849393957, 0.719718023452213, 0.678493086851215, 0.621074494722543,
    0.548389808139718, 0.461730945425636, 0.36273331202841, 0.253346732359107,
    0.135798653169166, 0.012550272501144, -0.113753577863184,
    -0.240339809358788, -0.36436658535247, -0.482985715635579,
    -0.593406597748258, -0.692960301259136, -0.779162356320678,
    -0.84977280460389, -0.902852099869886, -0.936811507032296,
    -0.950456741947632, -0.943023717954742, -0.91420541726803,
    -0.864169082927337, -0.793563126698113, -0.703513366103024,
    -0.595608435157227, -0.471874453474112, -0.334739281999959,
    -0.186986935322639, -0.031702954800321, 0.127788231766758,
    0.287992734816844, 0.445322002652177, 0.596171433954132, 0.737001120196412,
    0.864416949174725, 0.975250252948618, 1.06663417536608]
oscillator_euler_x2 = [
    1, 1, 0.9775, 0.9325, 0.86550625, 0.77753125, 0.670082359375,
    0.545139015625, 0.405118818789062, 0.252832994101562, 0.091431995991309,
    -0.075657744486231, -0.244804704873574, -0.412249366009977,
    -0.574185921286725, -0.726846865828249, -0.866588627140821,
    -0.989976333972257, -1.09386579669303, -1.17548079189942,
    -1.23248380668022, -1.26303850364328, -1.26586231495604, -1.24026775993682,
    -1.18619130283109, -1.10420882112679, -0.995537035108782,
    -0.862020550615424, -0.706104482832118, -0.530792952659966,
    -0.33959407162409, -0.136452349153366, 0.074330239928901,
    0.288183006867118, 0.500363343406935, 0.706059562292242, 0.900497605950893,
    1.07904930945797, 1.23733981683115, 1.37135171474152, 1.4775234667732,
    1.55283980522319, 1.59491186567078, 1.60204503050086, 1.57329267835334,
    1.50849431301955, 1.40829686242281, 1.27415828978313, 1.10833303773893,
    0.913839224174621, 0.694407917261181]
