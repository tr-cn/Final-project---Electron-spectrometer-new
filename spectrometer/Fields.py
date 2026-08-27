# import numpy as np

# class Magenetic_field:
#     def __init__(self, B0_T, width_mm, depth_mm, k=0, fringe=1, sharp_edge=1):
#         self.B0_T = B0_T
#         self.k = k # Gradient along Z direction
#         self.fringe = fringe
        
        # המרה למטרים עבור החישובים הפיזיקליים (כיוון שהאינטגרטור עובד במטרים)
#         self.D_m = width_mm * 1e-3
#         self.L_m = depth_mm * 1e-3
        
#         if sharp_edge:
#             self.a0 = 0.3835
#             self.a1 = 2.388
#             self.a2 = -0.8171
#             self.a3 = 0.2 
#         else:
#             self.a0 = 0.531
#             self.a1 = 2.341
#             self.a2 = 0.7799
#             self.a3 = 0.110
            
#     def _poly_S(self, s):
#         """ חישוב הפולינום עבור מרחק מנורמל s """
#         return self.a0 + self.a1 * s + self.a2 * (s**2) + self.a3 * (s**3)
    
#     def _enge_factor(self, s_complex):
#         """ חישוב פונקציית Enge בטוחה למספרים מרוכבים (מונעת קריסת Overflow) """
#         S = self._poly_S(s_complex)
        
        # אם אנחנו עמוק בתוך המגנט (S שלילי מאוד) או רחוק בחוץ (S חיובי מאוד)
#         if np.real(S) > 100:
#             return 0.0 + 0.0j
#         elif np.real(S) < -100:
#             return 1.0 + 0.0j
            
#         return 1.0 / (1.0 + np.exp(S))
    
#     def _get_magnetic_field(self, R):
#         if self.fringe == 0:
            # מודל ללא פרינג': שדה קיים רק בין הכניסה ליציאה
#             y = R[1]
#             if 0 <= y <= self.L_m:
#                 return np.array([self.B0_T, 0, 0])
#             return np.array([0, 0, 0])
            
#         x = R[0]
#         y = R[1]
#         z = R[2]
        
#         y_complex = y + 1j * x
        
#         # 1. Enge של אזור הכניסה (סביב y=0)
        # נרמול: בפנים y חיובי לכן נכפיל במינוס כדי לקבל s שלילי (F->1)
#         s_in = -y_complex / self.D_m
#         F_in = self._enge_factor(s_in)
        
#         # 2. Enge של אזור היציאה (סביב y=L)
        # נרמול: בפנים y קטן מ-L לכן (y-L) שלילי (F->1)
#         s_out = (y_complex - self.L_m) / self.D_m
#         F_out = self._enge_factor(s_out)
        
        # הפונקציה המלאה היא מכפלת הקצוות
#         F = F_in * F_out
        
#         Bx = self.B0_T * np.real(F) * (1 - self.k * z)
#         By = self.B0_T * np.imag(F)
#         Bz = self.B0_T * self.k * x
        
#         return np.array([Bx, By, Bz])


# def plot_magnetic_quiver(ax, spec, width_mm, depth_mm, height_mm):
#     """
#     מצייר את קווי השדה המגנטי בעזרת חיצים (Quiver) על גבי הגרף התלת-ממדי.
#     """
    # רזולוציית הרשת
#     nx, ny, nz = 5, 25, 5 
    
    # יצירת רשת הקואורדינטות במטרים
#     X_width = np.linspace(-width_mm/2, width_mm/2, nx) * 1e-3
#     Y_depth = np.linspace(-depth_mm * 0.3, depth_mm * 1.3, ny) * 1e-3 
#     Z_height = np.linspace(-height_mm/2, height_mm/2, nz) * 1e-3
    
#     X, Y, Z = np.meshgrid(X_width, Y_depth, Z_height, indexing='ij')
    
#     Bx = np.zeros_like(X)
#     By = np.zeros_like(X)
#     Bz = np.zeros_like(X)
    
    # חישוב השדה המגנטי בכל נקודה ברשת
#     for i in range(nx):
#         for j in range(ny):
#             for k in range(nz):
#                 R_current_m = np.array([X[i,j,k], Y[i,j,k], Z[i,j,k]])
#                 B_vec = spec._get_magnetic_field(R_current_m)
                
#                 Bx[i,j,k] = B_vec[0]
#                 By[i,j,k] = B_vec[1]
#                 Bz[i,j,k] = B_vec[2]
                
#     X_mm = X * 1e3
#     Y_mm = Y * 1e3
#     Z_mm = Z * 1e3
    
    # ציור החיצים בגרף
#     ax.quiver(X_mm, Y_mm, Z_mm, 
#               Bx, By, Bz, 
#               length=10.0,       # אורך החץ הוגדל כדי לפצות על ביטול הנרמול
#               normalize=False,   # עכשיו החיצים באמת יתקצרו היכן שהשדה דועך!
#               color='cyan', 
#               alpha=0.6, 
#               arrow_length_ratio=0.3)









import numpy as np


class Magenetic_field:
    def __init__(self, B0_T, width_mm=12.5, depth_mm=50.8, k=0, fringe = 1, sharp_edge=1, yoke=1):
        self.B0_T = B0_T
        self.k = k # Gradient along Z direction
        self.fringe = fringe
        self.yoke = yoke
        
        self.width_m = width_mm*1e-3
        self.depth_m = depth_mm*1e-3
        if sharp_edge:
            self.a0 = 0.3835
            self.a1 = 2.388
            self.a2 = -0.8171
            self.a3 = 0.2 
        else:
            self.a0 = 0.531
            self.a1 = 2.341
            self.a2 = 0.7799
            self.a3 = 0.110
        
        
    
        
    def _enge_func(self, y_complex):
        S = self.a0 + self.a1 * y_complex + self.a2 * (y_complex**2) + self.a3 * (y_complex**3)
        return 1 / (1 + np.exp(S))
    
    def _get_magnetic_field(self,R):
        if self.fringe == 0:
            return np.array([self.B0_T,0,0])
        
        
        x = R[0]; y = R[1]; z = R[2];
        
        y_complex = y + 1j*x
        y_comlex_norm = -y_complex/self.width_m
        F_in = self._enge_func(y_comlex_norm)
        
        if self.yoke:
            F = F_in
        else:
            y_complex_norm_out = (y_complex - self.depth_m)/self.width_m
            F_out = self._enge_func(y_complex_norm_out)
            F = F_in*F_out
        
        Bx = self.B0_T * np.real(F) * (1 - self.k * z)
        By = self.B0_T * np.imag(F)
        Bz = self.B0_T * self.k * x
        
        return np.array([Bx,By,Bz])
        

class Electric_field:
    def __init__(self, E0_V0m, width_mm=12.5, depth_mm=50.8, k=0, fringe = 1, sharp_edge=1, yoke=1):
        self.E0_V0m = E0_V0m
        self.k = k # Gradient along Z direction
        self.fringe = fringe
        self.yoke = yoke
        
        self.width_m = width_mm*1e-3
        self.depth_m = depth_mm*1e-3
        if sharp_edge:
            self.a0 = 0.3835
            self.a1 = 2.388
            self.a2 = -0.8171
            self.a3 = 0.2 
        else:
            self.a0 = 0.531
            self.a1 = 2.341
            self.a2 = 0.7799
            self.a3 = 0.110
        
        
    
        
    def _enge_func(self, y_complex):
        S = self.a0 + self.a1 * y_complex + self.a2 * (y_complex**2) + self.a3 * (y_complex**3)
        return 1 / (1 + np.exp(S))
    
    def _get_magnetic_field(self,R):
        if self.fringe == 0:
            return np.array([self.B0_T,0,0])
        
        
        x = R[0]; y = R[1]; z = R[2];
        
        y_complex = y + 1j*x
        y_comlex_norm = -y_complex/self.width_m
        F_in = self._enge_func(y_comlex_norm)
        
        if self.yoke:
            F = F_in
        else:
            y_complex_norm_out = (y_complex - self.depth_m)/self.width_m
            F_out = self._enge_func(y_complex_norm_out)
            F = F_in*F_out
        
        Bx = self.B0_T * np.real(F) * (1 - self.k * z)
        By = self.B0_T * np.imag(F)
        Bz = self.B0_T * self.k * x
        
        return np.array([Bx,By,Bz])

        
        
        
def plot_magnetic_quiver(ax, mag_T, width_mm, depth_mm, height_mm, yoke,fringe, shield_mm):

    
    quive_len = 2;
# Gread resolution
    nx, ny, nz = 5, 25, 5 
    
# Gread generation
    X_width = np.linspace(-width_mm/2, width_mm/2 - quive_len, nx) * 1e-3
    if yoke:
        if fringe:
            Y_depth = np.linspace(-shield_mm, depth_mm * 1, ny) * 1e-3 
        else:
            Y_depth = np.linspace(0, depth_mm * 1, ny) * 1e-3 
    else:
        if fringe:
            Y_depth = np.linspace(-shield_mm, depth_mm * 1.3, ny) * 1e-3 
        else:
            Y_depth = np.linspace(0, depth_mm, ny) * 1e-3 
            
    Z_height = np.linspace(-height_mm/2, height_mm/2, nz) * 1e-3
    
    X, Y, Z = np.meshgrid(X_width, Y_depth, Z_height, indexing='ij')
    
    Bx = np.zeros_like(X)
    By = np.zeros_like(X)
    Bz = np.zeros_like(X)
    

# Caclulating the magnetic field in every point
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                R_current_m = np.array([X[i,j,k], Y[i,j,k], Z[i,j,k]])
                B_vec = mag_T._get_magnetic_field(R_current_m)
                
                Bx[i,j,k] = B_vec[0]
                By[i,j,k] = B_vec[1]
                Bz[i,j,k] = B_vec[2]
                
    X_mm = X * 1e3
    Y_mm = Y * 1e3
    Z_mm = Z * 1e3
    
# plot magentic filed arrows

    ax.quiver(X_mm, Y_mm, Z_mm, 
              Bx, By, Bz, 
              length=quive_len,       
              normalize=False,
              color='cyan', 
              alpha=0.6, 
              arrow_length_ratio=0.3)