# Newton_48.py
import math

def newton_system_48(x0, y0, tol=1e-6, max_iter=50):
    x, y = x0, y0
    print("=== Newton-Raphson (sistem) ===")
    print(f"{'i':<4}{'x':<14}{'y':<14}{'Δx':<12}{'Δy':<12}")
    print(f"{0:<4}{x:<14.6f}{y:<14.6f}{0.0:<12.6f}{0.0:<12.6f}")
    for i in range(1, max_iter+1):
        f1 = x*x + x*y - 10
        f2 = y + 3*x*(y**2) - 57
        J11 = 2*x + y
        J12 = x
        J21 = 3*(y**2)
        J22 = 1 + 6*x*y
        det = J11*J22 - J12*J21
        if abs(det) < 1e-14:
            print("Jacobian near singular")
            return x, y
        dx = - ( f1*J22 - f2*J12 ) / det
        dy = - ( f2*J11 - f1*J21 ) / det
        x_new = x + dx
        y_new = y + dy
        print(f"{i:<4}{x_new:<14.6f}{y_new:<14.6f}{abs(dx):<12.6f}{abs(dy):<12.6f}")
        if abs(dx) < tol and abs(dy) < tol:
            print("Konvergen pada iterasi ke-", i)
            return x_new, y_new
        x, y = x_new, y_new
    print("Tidak konvergen.")
    return x, y

# run example
if __name__ == "__main__":
    newton_system_48(1.5, 3.5)
