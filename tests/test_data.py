'''Defines fuctions and data for tests.'''

__all__ = [
        'oscillator_1st_deriv',
        'oscillator_2nd_deriv',
        'oscillator_euler_t',
        'oscillator_euler_x1',
        'oscillator_euler_x2',
        'oscillator_verlet_t',
        'oscillator_verlet_x1',
        'oscillator_verlet_v1',
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
    x = X[0]
    ddx = -(w**2) * x
    ddX = [ddx]
    return ddX

# Test data generated in a spreadsheet

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

oscillator_verlet_t = oscillator_euler_t
oscillator_verlet_x1 = [
    0, 0.1, 0.19775, 0.291050625, 0.3778026109375, 0.456054038128906,
    0.524044249462412, 0.580243465183014, 0.623387202936998, 0.652504728624899,
    0.66694089791874, 0.66637089700941, 0.650807550917367, 0.620601034929685,
    0.576430995656084, 0.519291258980221, 0.450467468977304, 0.371508160922397,
    0.284189919246736, 0.190477404388024, 0.092479147930581, -0.0075998893553,
    -0.107507929130687, -0.204997040500633, -0.297873718459315,
    -0.384048237752662, -0.461581671696575, -0.528729518027314,
    -0.58398095020244, -0.62609281099801, -0.654117583546125,
    -0.667424710464452, -0.665714781397329, -0.649026269748766,
    -0.617734667030856, -0.572544034304751, -0.51447116080679,
    -0.444822686190676, -0.365165701135271, -0.277292487804324,
    -0.183180193497778, -0.084946344837533, 0.015198796581557,
    0.115001965077561, 0.212217589359321, 0.304658317880495, 0.390244234249359,
    0.467049655347612, 0.533346459200544, 0.587642967721463, 0.62871750946865]
oscillator_verlet_v1 = [
    1, 0.98875, 0.955253125, 0.9002630546875, 0.825017065644531,
    0.731208192624561, 0.620947135270537, 0.496714767372927, 0.361306317209426,
    0.217768474908713, 0.069330841922554, -0.080666735006863,
    -0.228849310398625, -0.371882776306418, -0.506548879747317,
    -0.629817633393901, -0.738915490289122, -0.831387748652838,
    -0.905153782671865, -0.958553856580776, -0.990386468716619,
    -0.999935385306338, -0.986985755726665, -0.951828946643141,
    -0.895255986260147, -0.8185397661863, -0.723406401373261,
    -0.611996392529323, -0.486816464853476, -0.350683166718426,
    -0.206659497332211, -0.05798598925602, 0.09199220357843, 0.239900571832366,
    0.382411177220073, 0.516317531120329, 0.638606740570377, 0.746527298357592,
    0.837650991931761, 0.909927538187465, 0.961730714833952, 0.991894950396674,
    0.999741549575472, 0.985093963888821, 0.948281764014672, 0.890133224450192,
    0.811956687335583, 0.715511124755924, 0.602966561869257, 0.476855251340532,
    0.340014697656644]