import numpy as np
from math import *
import graphics as gr

g = 9.8
a = 0.25
v = 1
l = 2
h = 0.01
start_acceleration = 1.5

pixel_per_meter = 150
length = l*pixel_per_meter

def calculateF(t, y):
    return -(g+a*v*v*cos(v*t))*sin(y)/l

def update_coords(coords, angle, center_point):
    new_point = gr.Point(center_point.getX() + length * sin(angle),
                         center_point.getY() + length * cos(angle))
    velocity = gr.Point(new_point.getX() - coords.getX(),
                        new_point.getY() - coords.getY())
    return velocity

def paint(list, method):
    SIZE_X = 800
    SIZE_Y = 800
    window = gr.GraphWin("Маятник Капіци ({})".format(method), SIZE_X, SIZE_Y)
    angle = 0 + pi
    center_point = gr.Point(SIZE_X / 2, SIZE_Y / 2)
    coords = gr.Point(center_point.getX() + length * sin(angle),
                      center_point.getY() + length * cos(angle))
    # Фон
    rectangle = gr.Rectangle(gr.Point(0, 0), gr.Point(SIZE_X, SIZE_Y))
    rectangle.setFill('white')
    rectangle.draw(window)
    # центр
    center = gr.Circle(center_point, 7)
    center.setFill('black')
    center.draw(window)
    # куля
    ball = gr.Circle(gr.Point(coords.getX(), coords.getY()), 25)
    ball.setFill('green')
    ball.draw(window)

    ceil_point = gr.Point(SIZE_X / 2, SIZE_Y / 2 - a * pixel_per_meter)
    float_point = gr.Point(SIZE_X / 2, SIZE_Y / 2 + a * pixel_per_meter)

    line_pidvis = gr.Line(ceil_point, float_point)
    line_pidvis.setFill('black')
    line_pidvis.draw(window)

    sign = True ## в яку сторону рухатися центру
    temp = 0

    new_y = 0
    diff_y = 3 * a * v / 100

    for item in list:
        old_value = SIZE_Y / 2 + new_y * pixel_per_meter
        ## поміняти центр
        if sign == True:
            new_y = new_y + diff_y + temp
            if new_y >= a:
                sign = False
                center_point = gr.Point(center_point.getX(), SIZE_Y / 2 + a * pixel_per_meter)
                center.move(0, -center_point.getY() + old_value)
                temp = new_y - a
                new_y = a
            else:
                center_point = gr.Point(center_point.getX(), SIZE_Y / 2 + new_y * pixel_per_meter)
                center.move(0, -center_point.getY() + old_value)
                temp = 0
        else:
            new_y = new_y - diff_y - temp
            if new_y <= -a:
                sign = True
                center_point = gr.Point(center_point.getX(), SIZE_Y / 2 - a * pixel_per_meter)
                center.move(0, -center_point.getY() + old_value)
                temp = -a - new_y
                new_y = -a
            else:
                center_point = gr.Point(center_point.getX(), SIZE_Y / 2 + new_y * pixel_per_meter)
                center.move(0, -center_point.getY() + old_value)
                temp = 0
        ###
        if not window.isClosed():
            angle = item + pi
            velocity = update_coords(ball.getCenter(), angle, center.getCenter())
            ball.move(velocity.x, velocity.y)
            line = gr.Line(center.getCenter(), ball.getCenter())
            line.draw(window)
            gr.time.sleep(0.01)
            line.undraw()


    line_pidvis.undraw()
    window.close()

def calculate_RK_Adams():
    y_curr = 0
    z_curr = start_acceleration

    resZ_rk = []
    res_rk = []

    resZ_rk.append(z_curr)
    res_rk.append(y_curr)

    TIME = 20

    for x_curr in np.arange(0, TIME, h):
        k0 = z_curr
        q0 = calculateF(x_curr, y_curr)

        k1 = z_curr + q0 * h / 2
        q1 = calculateF(x_curr + h / 2, y_curr + k0 * h / 2)

        k2 = z_curr + q1 * h / 2
        q2 = calculateF(x_curr + h / 2, y_curr + k1 * h / 2)

        k3 = z_curr + q2 * h
        q3 = calculateF(x_curr + h, y_curr + k2 * h)

        y_next = y_curr + h / 6 * (k0 + 2 * (k1 + k2) + k3)
        z_next = z_curr + h / 6 * (q0 + 2 * (q1 + q2) + q3)

        resZ_rk.append(z_next)
        res_rk.append(y_next)

        y_curr = y_next
        z_curr = z_next

    # Adams
    x_prev = 0
    x_curr = 0.5
    x_next = 1

    y_prev = 0.05
    z_prev = 1.5

    z_curr = resZ_rk[1]
    y_curr = res_rk[1]

    z_next = resZ_rk[2]
    y_next = res_rk[2]

    res_adams = []
    res_adams.append(y_prev)
    res_adams.append(y_curr)
    res_adams.append(y_next)

    time = 3
    while time < TIME / h + 1:
        v_prev = calculateF(x_prev, y_prev)
        v_curr = calculateF(x_curr, y_curr)
        v_next = calculateF(x_next, y_next)

        deltaZ = h * v_next + 0.5 * h * (v_next - v_curr) + 5 / 12 * h * ((v_next - v_curr) - (v_curr - v_prev))
        z_res = z_next + deltaZ  # Zn+1
        deltaY = h * z_next + 0.5 * h * (z_next - z_curr) + 5 / 12 * h * ((z_next - z_curr) - (z_curr - z_prev))
        y_res = y_next + deltaY
        res_adams.append(y_res)

        y_prev = y_next
        y_curr = y_next
        y_next = y_res

        z_prev = z_next
        z_curr = z_next
        z_next = z_res

        x_prev += h
        x_curr += h
        time += 1
    return res_rk, res_adams

if __name__ == "__main__":
    res_rk, res_adams = calculate_RK_Adams()
    while(True):
        print("Choose the method: 1.Runge-Kutt 2. Adams")
        choose = input()
        if choose == "1":
            paint(res_rk, "RK")
        elif choose == "2":
            paint(res_adams, "Adams")
        else:
            break
