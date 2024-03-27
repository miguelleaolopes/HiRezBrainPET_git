from calendar import c
import re
import tkinter as tk
from tkinter import ttk
from wsgiref import validate
import ttkbootstrap as ttk
from tkinter import messagebox
from tkinter import filedialog
from ttkbootstrap import Style
import os
import numpy as np
import locale

# Set the locale to C at the start of your script
locale.setlocale(locale.LC_NUMERIC, 'C')

class BatchScriptGenerator(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Python GUI for Running CASToR and Script Generator")
        self.script_dir = os.path.dirname(os.path.realpath(__file__))
        # Png logo
        self.iconphoto(False, tk.PhotoImage(file=os.path.join(self.script_dir, 'CASToR_logo.png')))
        self.resizable(True, True)
        self.minsize(675, 482)
        self.style = Style(theme='darkly')

        # Variables to store user inputs
        self.mpi_bool_var = tk.BooleanVar()
        self.mpi_threads_var = tk.IntVar()
        self.verbose_level_var = tk.IntVar()
        self.last_iter_bool_var = tk.BooleanVar()
        self.flip_var = tk.StringVar()
        self.stats_need_bool_var = tk.BooleanVar()

        self.main_program_path_var = tk.StringVar()
        self.datafile_path_var = tk.StringVar()
        self.output_path_var = tk.StringVar()
        self.sensitivity_path_var = tk.StringVar()
        self.configuration_path_var = tk.StringVar()

        self.voxel_number_var = tk.StringVar()
        self.voxel_size_var = tk.StringVar()
        self.fov_size_var = tk.StringVar()
        self.offset_var = tk.StringVar()
        self.voxel_size_x_var = tk.DoubleVar()
        self.voxel_size_y_var = tk.DoubleVar()
        self.voxel_size_z_var = tk.DoubleVar()
        self.voxel_number_x_var = tk.IntVar()
        self.voxel_number_y_var = tk.IntVar()
        self.voxel_number_z_var = tk.IntVar()
        self.fov_size_x_var = tk.DoubleVar()
        self.fov_size_y_var = tk.DoubleVar()
        self.fov_size_z_var = tk.DoubleVar()
        self.offset_x_var = tk.DoubleVar()
        self.offset_y_var = tk.DoubleVar()
        self.offset_z_var = tk.DoubleVar()

        self.iterations_var = tk.StringVar()
        self.optimizer_var = tk.StringVar()
        self.projector_var = tk.StringVar()
        self.penalty_var = tk.StringVar()
        self.penalty_strength_var = tk.DoubleVar()
        self.convolution_need_bool_var = tk.BooleanVar()
        self.convolution_num_var = tk.IntVar()
        self.previous_num_conv = int()
        self.convolution_type_vars = []
        self.convolution_value_vars = []
        self.convolution_x_var = []
        self.convolution_y_var = []
        self.convolution_sigma_var = []


        # Set initial values
        self.set_initial_values()

        # Create GUI elements
        self.create_widgets()

    def set_initial_values(self):
        # Set initial default values for the parameters
        self.mpi_bool_var.set(True)
        self.mpi_threads_var.set(0)
        self.verbose_level_var.set(2)
        self.last_iter_bool_var.set(False)
        self.flip_var.set("None")
        self.stats_need_bool_var.set(True)

        self.main_program_path_var.set(r"..\castor-recon.exe")
        self.datafile_path_var.set(r"..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh")
        self.output_path_var.set(r"..\..\Results")
        self.sensitivity_path_var.set("")
        self.configuration_path_var.set("")

        self.voxel_number_var.set("0.2,0.2,0.4")
        self.voxel_size_var.set("200,200,50")
        self.fov_size_var.set("40,40,20")
        self.offset_var.set("0,-45,-8")
        self.voxel_size_x_var.set(0.2)
        self.voxel_size_y_var.set(0.2)
        self.voxel_size_z_var.set(0.4)
        self.voxel_number_x_var.set(200)
        self.voxel_number_y_var.set(200)
        self.voxel_number_z_var.set(50)
        self.fov_size_x_var.set(40)
        self.fov_size_y_var.set(40)
        self.fov_size_z_var.set(20)
        self.offset_x_var.set(0)
        self.offset_y_var.set(-45)
        self.offset_z_var.set(-8)

        self.iterations_var.set("2:50,2:10,5:5")
        self.optimizer_var.set("MLEM")
        self.projector_var.set("joseph")
        self.penalty_var.set("MRF")
        self.penalty_strength_var.set(0.5)
        self.convolution_need_bool_var.set(True)
        self.convolution_num_var.set(1)
        self.convolution_type_vars = []
        self.convolution_value_vars = []
        self.convolution_x_var = []
        self.convolution_y_var = []
        self.convolution_sigma_var = []
        for i in range(self.convolution_num_var.get()):
            type_var = tk.StringVar()
            type_var.set("psf")  # Set initial value
            value_var = tk.StringVar()
            value_var.set("gaussian,1,1,3::psf")  # Set initial value
            value_x_var = tk.DoubleVar()
            value_x_var.set(1.0)  # Set initial value
            value_y_var = tk.DoubleVar()
            value_y_var.set(1.0)  # Set initial value
            value_sigma_var = tk.DoubleVar()
            value_sigma_var.set(3.0)  # Set initial value
            self.convolution_type_vars.append(type_var)
            self.convolution_value_vars.append(value_var)
            self.convolution_x_var.append(value_x_var)
            self.convolution_y_var.append(value_y_var)
            self.convolution_sigma_var.append(value_sigma_var)

    def update_entries(self):
        # Update the entries that depend on other widgets
        self.voxel_number_var.set(f"{self.voxel_size_x_var.get()},{self.voxel_size_y_var.get()},{self.voxel_size_z_var.get()}")
        self.voxel_size_var.set(f"{self.voxel_number_x_var.get()},{self.voxel_number_y_var.get()},{self.voxel_number_z_var.get()}")
        self.fov_size_var.set(f"{self.fov_size_x_var.get()},{self.fov_size_y_var.get()},{self.fov_size_z_var.get()}")
        self.offset_var.set(f"{self.offset_x_var.get()},{self.offset_y_var.get()},{self.offset_z_var.get()}")
        self.convolution_value_vars = [f"gaussian,{self.convolution_x_var[i].get()},{self.convolution_y_var[i].get()},{self.convolution_sigma_var[i].get()}::{self.convolution_type_vars[i].get()}" for i in range(self.convolution_num_var.get())]
        # NEED TO CREATE self.button etc...

    def create_widgets(self):
        # Add GUI elements for user input
        options_frame = ttk.Frame(self)
        options_frame.pack(padx=10, pady=10)
        ttk.Sizegrip(self, style='info.TSizegrip').pack(side="right", fill="y")

        # Add a menubar
        self.create_menu()

        # Path variables section
        self.create_path_variables_widgets(options_frame)

        # Image specifications section
        self.create_image_specifications_widgets(options_frame)

        # Miscellaneous section
        self.create_miscellaneous_widgets(options_frame)

        # Algorithm section
        self.create_algorithm_widgets(options_frame)

        # Convolution section
        self.create_convolution_widgets(options_frame)

        button_frame = ttk.Frame(self)
        button_frame.pack(padx=10, pady=(0, 15))

        # Button to generate batch script
        generate_button = ttk.Button(button_frame, text="Generate batch script", command=self.generate_script)
        generate_button.pack(side="left", padx=5)

        # Button to run CASToR
        run_button = ttk.Button(button_frame, text="Run CASToR program", command=self.run_castor_program)
        run_button.pack(side="left", padx=5)
    
    def create_menu(self):
        menubar = tk.Menu(self)
        self.config(menu=menubar)
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Reset variable to default", command=lambda: (self.set_initial_values(),self.update_convolution_entries()))
        file_menu.add_command(label="Print all variables", command=self.print_test_all_variables)
        file_menu.add_command(label="Generate Script", command=self.generate_script)
        file_menu.add_command(label="Run CASToR Program", command=self.run_castor_program)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.quit)
        # Add a open from file menu
        menubar.add_command(label="Open from file", command=self.open_from_file)
        # Add a aboutme menu
        menubar.add_command(label="About Me", command=self.show_about)

    def create_path_variables_widgets(self, options_frame):
        path_frame = ttk.LabelFrame(options_frame, text="Path Variables")
        path_frame.grid(row=0, column=0, padx=10, pady=5, sticky="nsew")

        self.folder_icon = tk.PhotoImage(file=os.path.join(self.script_dir, 'folder.png'))
        self.folder_icon = self.folder_icon.subsample(8, 8)

        ttk.Label(path_frame, text="CASToR Main Prog:").grid(row=0, column=0, sticky="nsew", padx=(5, 2))
        ttk.Entry(path_frame, textvariable=self.main_program_path_var,validate="focusout",validatecommand=lambda: self.validate_program_path(self.main_program_path_var)).grid(row=0, column=1, pady=3)
        ttk.Button(path_frame, image=self.folder_icon, command=lambda: self.browse_path(self.main_program_path_var)).grid(row=0, column=2, pady=3, padx=(0, 5))
        
        ttk.Label(path_frame, text="Datafile Cdh Path:").grid(row=1, column=0, sticky="nsew", padx=(5, 2))
        ttk.Entry(path_frame, textvariable=self.datafile_path_var,validate="focusout",validatecommand=lambda: self.validate_program_path(self.datafile_path_var)).grid(row=1, column=1, pady=3)
        ttk.Button(path_frame, image=self.folder_icon, command=lambda: self.browse_path(self.datafile_path_var)).grid(row=1, column=2, pady=3, padx=(0, 5))

        ttk.Label(path_frame, text="Output Folder:").grid(row=2, column=0, sticky="nsew", padx=(5, 2))
        ttk.Entry(path_frame, textvariable=self.output_path_var).grid(row=2, column=1, pady=3)
        ttk.Button(path_frame, image=self.folder_icon, command=lambda: self.browse_path(self.output_path_var)).grid(row=2, column=2, pady=3, padx=(0, 5))

        ttk.Label(path_frame, text="Sensitivity Path:").grid(row=3, column=0, sticky="nsew", padx=(5, 2))
        ttk.Entry(path_frame, textvariable=self.sensitivity_path_var,validate="focusout",validatecommand=lambda: self.validate_program_path(self.sensitivity_path_var)).grid(row=3, column=1, pady=(3, 6))
        ttk.Button(path_frame, image=self.folder_icon, command=lambda: self.browse_path(self.sensitivity_path_var)).grid(row=3, column=2, pady=(3, 6), padx=(0, 5))

    def create_image_specifications_widgets(self, options_frame):
        image_frame = ttk.LabelFrame(options_frame, text="Image Specifications (X,Y,Z)")
        image_frame.grid(row=0, column=1, padx=10, pady=5, sticky="nsew")

        def set_voxel_size(value, axis):
            if axis == 'x':
                self.voxel_size_x_var.set(value)
            elif axis == 'y':
                self.voxel_size_y_var.set(value)
            elif axis == 'z':
                self.voxel_size_z_var.set(value)
            update_fov()

        def set_voxel_number(value, axis):
            if axis == 'x':
                self.voxel_number_x_var.set(value)
            elif axis == 'y':
                self.voxel_number_y_var.set(value)
            elif axis == 'z':
                self.voxel_number_z_var.set(value)
            update_fov()

        def set_fov_size(value, axis):
            if axis == 'x':
                self.fov_size_x_var.set(value)
            elif axis == 'y':
                self.fov_size_y_var.set(value)
            elif axis == 'z':
                self.fov_size_z_var.set(value)
            update_voxel_size()

        def set_offset(value, axis):
            if axis == 'x':
                self.offset_x_var.set(value)
            elif axis == 'y':
                self.offset_y_var.set(value)
            elif axis == 'z':
                self.offset_z_var.set(value)

        def validate_voxel_size(P, axis):
            if P == "":
                value = 0.2
            elif P.isdigit() or self.is_valid_float(P):
                value = round(float(P),5)
                if value > 1000:
                    value = 1000
                elif value <= 0:
                    value = 0.2
            else:
                value = 0.2

            self.after(10, lambda: set_voxel_size(value, axis))
            return True

        def validate_voxel_number(P, axis):
            if P == "":
                value = 200
            elif P.isdigit() or self.is_valid_float(P):
                value = round(float(P))
                if value > 100000:
                    value = 100000
                elif value <= 0:
                    value = 1
            else:
                value = 200

            self.after(10, lambda: set_voxel_number(value, axis))
            return True

        def validate_fov_size(P, axis):
            # P is the new value of the Spinbox after the edit
            if P == "":
                value = 40
            elif P.isdigit() or self.is_valid_float(P):
                value = round(float(P),5)
                if value > 1000000000:
                    value = 1000000000
                elif value <= 0:
                    value = 1
            else:
                value = 40

            self.after(10, lambda: set_fov_size(value, axis))
            return True

        def validate_offset(P, axis):
            # P is the new value of the Spinbox after the edit
            if P == "":
                value = 0
            elif P.isdigit() or self.is_valid_float(P):
                value = round(float(P),5)
                if value > 1000000:
                    value = 1000000
                elif value < -1000000:
                    value = -1000000
            else:
                value = 0

            self.after(10, lambda: set_offset(value, axis))
            return True
        
        vcmd_voxel_number_x = (self.register(lambda P: validate_voxel_number(P, 'x')), '%P')
        vcmd_voxel_number_y = (self.register(lambda P: validate_voxel_number(P, 'y')), '%P')
        vcmd_voxel_number_z = (self.register(lambda P: validate_voxel_number(P, 'z')), '%P')
        vcmd_voxel_size_x = (self.register(lambda P: validate_voxel_size(P, 'x')), '%P')
        vcmd_voxel_size_y = (self.register(lambda P: validate_voxel_size(P, 'y')), '%P')
        vcmd_voxel_size_z = (self.register(lambda P: validate_voxel_size(P, 'z')), '%P')
        vcmd_fov_size_x = (self.register(lambda P: validate_fov_size(P, 'x')), '%P')
        vcmd_fov_size_y = (self.register(lambda P: validate_fov_size(P, 'y')), '%P')
        vcmd_fov_size_z = (self.register(lambda P: validate_fov_size(P, 'z')), '%P')
        vcmd_offset_x = (self.register(lambda P: validate_offset(P, 'x')), '%P')
        vcmd_offset_y = (self.register(lambda P: validate_offset(P, 'y')), '%P')
        vcmd_offset_z = (self.register(lambda P: validate_offset(P, 'z')), '%P')

        def update_fov(*args):
            voxel_size_x = self.voxel_size_x_var.get() if self.voxel_size_x_var.get() else 0.
            voxel_size_y = self.voxel_size_y_var.get() if self.voxel_size_y_var.get() else 0
            voxel_size_z = self.voxel_size_z_var.get() if self.voxel_size_z_var.get() else 0
            voxel_number_x = self.voxel_number_x_var.get() if self.voxel_number_x_var.get() else 0
            voxel_number_y = self.voxel_number_y_var.get() if self.voxel_number_y_var.get() else 0
            voxel_number_z = self.voxel_number_z_var.get() if self.voxel_number_z_var.get() else 0
            fov_size_x = voxel_size_x * voxel_number_x
            fov_size_y = voxel_size_y * voxel_number_y
            fov_size_z = voxel_size_z * voxel_number_z
            self.fov_size_x_var.set(round(fov_size_x,5))
            self.fov_size_y_var.set(round(fov_size_y,5))
            self.fov_size_z_var.set(round(fov_size_z,5))

        def update_voxel_size(*args):
            fov_size_x = self.fov_size_x_var.get() if self.fov_size_x_var.get() else 0
            fov_size_y = self.fov_size_y_var.get() if self.fov_size_y_var.get() else 0
            fov_size_z = self.fov_size_z_var.get() if self.fov_size_z_var.get() else 0
            voxel_number_x = self.voxel_number_x_var.get() if self.voxel_size_x_var.get() else 0
            voxel_number_y = self.voxel_number_y_var.get() if self.voxel_size_y_var.get() else 0
            voxel_number_z = self.voxel_number_z_var.get() if self.voxel_size_z_var.get() else 0
            voxel_size_x = fov_size_x / voxel_number_x if voxel_number_x != 0 else 0
            voxel_size_y = fov_size_y / voxel_number_y if voxel_number_y != 0 else 0
            voxel_size_z = fov_size_z / voxel_number_z if voxel_number_z != 0 else 0
            self.voxel_size_x_var.set(round(voxel_size_x,5))
            self.voxel_size_y_var.set(round(voxel_size_y,5))
            self.voxel_size_z_var.set(round(voxel_size_z,5))

        ttk.Label(image_frame, text="Voxel Size (mm):").grid(row=0, column=0, sticky="nsew", padx=(5, 2))
        spinbox_voxsize_x = ttk.Spinbox(image_frame, textvariable=self.voxel_size_x_var, from_=0.0, to=1000.0,increment=0.1, width=6, validate='focusout', validatecommand=vcmd_voxel_size_x, command=update_fov)
        spinbox_voxsize_x.grid(row=0, column=1, pady=3)
        spinbox_voxsize_x.bind('<Return>', lambda event: validate_voxel_size(spinbox_voxsize_x.get(), 'x'))
        spinbox_voxsize_y = ttk.Spinbox(image_frame, textvariable=self.voxel_size_y_var, from_=0, to=1000,increment=0.1, width=6, validate='focusout', validatecommand=vcmd_voxel_size_y, command=update_fov)
        spinbox_voxsize_y.grid(row=0, column=2, pady=3)
        spinbox_voxsize_y.bind('<Return>', lambda event: validate_voxel_size(spinbox_voxsize_y.get(), 'y'))
        spinbox_voxsize_z = ttk.Spinbox(image_frame, textvariable=self.voxel_size_z_var, from_=0, to=1000,increment=0.1, width=6, validate='focusout', validatecommand=vcmd_voxel_size_z, command=update_fov)
        spinbox_voxsize_z.grid(row=0, column=3, pady=3, padx=(0, 5))
        spinbox_voxsize_z.bind('<Return>', lambda event: validate_voxel_size(spinbox_voxsize_z.get(), 'z'))

        ttk.Label(image_frame, text="Voxel Number:").grid(row=1, column=0, sticky="nsew", padx=(5, 2))
        spinbox_voxnb_x = ttk.Spinbox(image_frame, textvariable=self.voxel_number_x_var, from_=0, to=100000,increment=1, width=6, validate='focusout', validatecommand=vcmd_voxel_number_x, command=update_fov)
        spinbox_voxnb_x.grid(row=1, column=1, pady=3)
        spinbox_voxnb_x.bind('<Return>', lambda event: validate_voxel_number(spinbox_voxnb_x.get(), 'x'))
        spinbox_voxnb_y = ttk.Spinbox(image_frame, textvariable=self.voxel_number_y_var, from_=0, to=100000,increment=1, width=6, validate='focusout', validatecommand=vcmd_voxel_number_y, command=update_fov)
        spinbox_voxnb_y.grid(row=1, column=2, pady=3)
        spinbox_voxnb_y.bind('<Return>', lambda event: validate_voxel_number(spinbox_voxnb_y.get(), 'y'))
        spinbox_voxnb_z = ttk.Spinbox(image_frame, textvariable=self.voxel_number_z_var, from_=0, to=100000,increment=1, width=6, validate='focusout', validatecommand=vcmd_voxel_number_z, command=update_fov)
        spinbox_voxnb_z.grid(row=1, column=3, pady=3, padx=(0, 5))
        spinbox_voxnb_z.bind('<Return>', lambda event: validate_voxel_number(spinbox_voxnb_z.get(), 'z'))

        ttk.Label(image_frame, text="FOV Size (mm):").grid(row=2, column=0, sticky="nsew", padx=(5, 2))
        spinbox_fovsize_x = ttk.Spinbox(image_frame, textvariable=self.fov_size_x_var, from_=0, to=1000000000,increment=1, width=6, validate='focusout', validatecommand=vcmd_fov_size_x, command=update_voxel_size)
        spinbox_fovsize_x.grid(row=2, column=1, pady=3)
        spinbox_fovsize_x.bind('<Return>', lambda event: validate_fov_size(spinbox_fovsize_x.get(), 'x'))
        spinbox_fovsize_y = ttk.Spinbox(image_frame, textvariable=self.fov_size_y_var, from_=0, to=1000000000,increment=1, width=6, validate='focusout', validatecommand=vcmd_fov_size_y, command=update_voxel_size)
        spinbox_fovsize_y.grid(row=2, column=2, pady=3)
        spinbox_fovsize_y.bind('<Return>', lambda event: validate_fov_size(spinbox_fovsize_y.get(), 'y'))
        spinbox_fovsize_z = ttk.Spinbox(image_frame, textvariable=self.fov_size_z_var, from_=0, to=1000000000,increment=1, width=6, validate='focusout', validatecommand=vcmd_fov_size_z, command=update_voxel_size)
        spinbox_fovsize_z.grid(row=2, column=3, pady=3, padx=(0, 5))
        spinbox_fovsize_z.bind('<Return>', lambda event: validate_fov_size(spinbox_fovsize_z.get(), 'z'))

        ttk.Label(image_frame, text="Offset (mm):").grid(row=3, column=0, sticky="nsew", padx=(5, 2))
        spinbox_offset_x = ttk.Spinbox(image_frame, textvariable=self.offset_x_var, from_=-1000000,to=1000000, increment=1, width=6, validate='focusout', validatecommand=vcmd_offset_x)
        spinbox_offset_x.grid(row=3, column=1, pady=(3, 6))
        spinbox_offset_x.bind('<Return>', lambda event: validate_offset(spinbox_offset_x.get(), 'x'))
        spinbox_offset_y = ttk.Spinbox(image_frame, textvariable=self.offset_y_var, from_=-1000000,to=1000000, increment=1, width=6, validate='focusout', validatecommand=vcmd_offset_y)
        spinbox_offset_y.grid(row=3, column=2, pady=(3, 6))
        spinbox_offset_y.bind('<Return>', lambda event: validate_offset(spinbox_offset_y.get(), 'y'))
        spinbox_offset_z = ttk.Spinbox(image_frame, textvariable=self.offset_z_var, from_=-1000000,to=1000000, increment=1, width=6, validate='focusout', validatecommand=vcmd_offset_z)
        spinbox_offset_z.grid(row=3, column=3, pady=(3, 6), padx=(0, 5))
        spinbox_offset_z.bind('<Return>', lambda event: validate_offset(spinbox_offset_z.get(), 'z'))

    def create_miscellaneous_widgets(self, options_frame):
        misc_frame = ttk.LabelFrame(options_frame, text="Miscellaneous")
        misc_frame.grid(row=1, column=0, padx=10, pady=5, sticky="nsew")

        def validate_spinbox_mpi_input(P):
            # P is the new value of the Spinbox after the edit
            if P == "":
                self.mpi_threads_var.set(0)
                return True
            elif P.isdigit():
                if int(P) > 1000:
                    self.mpi_threads_var.set(1000)
                    return True
                elif int(P) < 0:
                    self.mpi_threads_var.set(0)
                    return True
                else:
                    return True
            else:
                self.mpi_threads_var.set(0)
                return True

        def validate_spinbox_verbose_input(P):
            # P is the new value of the Spinbox after the edit
            if P == "":
                self.verbose_level_var.set(2)
                return True
            elif P.isdigit():
                if int(P) > 5:
                    self.verbose_level_var.set(5)
                    return True
                elif int(P) < 1:
                    self.verbose_level_var.set(1)
                    return True
                else:
                    return True
            else:
                self.verbose_level_var.set(2)
                return True
        vcmd_mpi = (self.register(validate_spinbox_mpi_input), '%P')
        vcmd_verbose = (self.register(validate_spinbox_verbose_input), '%P')

        
        # Checkbuttons
        ttk.Checkbutton(misc_frame, text="MPI", variable=self.mpi_bool_var, style='Roundtoggle.Toolbutton').grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        ttk.Checkbutton(misc_frame, text="Print Stats", variable=self.stats_need_bool_var,style='Roundtoggle.Toolbutton').grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        ttk.Checkbutton(misc_frame, text="Save Last It", variable=self.last_iter_bool_var,style='Roundtoggle.Toolbutton').grid(row=2, column=0, sticky="nsew", padx=5, pady=5)
        # Add a separator line
        separator = ttk.Separator(misc_frame, orient="vertical")
        separator.grid(row=0, column=1, rowspan=3,sticky="ns", padx=10, pady=5)
        # Spinboxes and Menubutton
        ttk.Label(misc_frame, text="Nb of Threads:").grid(row=0, column=2, sticky="nsew")
        mpi_spinbox = ttk.Spinbox(misc_frame, textvariable=self.mpi_threads_var, from_=0,to=1000, increment=1, width=6, validate='focusout', validatecommand=vcmd_mpi)
        mpi_spinbox.grid(row=0, column=3)
        mpi_spinbox.bind('<Return>', lambda event: validate_spinbox_mpi_input(mpi_spinbox.get()))
        ttk.Label(misc_frame, text="Verbose Level:").grid(row=1, column=2, sticky="nsew")
        verbose_spinbox = ttk.Spinbox(misc_frame, textvariable=self.verbose_level_var, from_=1,to=5, increment=1, width=6, validate='focusout', validatecommand=vcmd_verbose)
        verbose_spinbox.grid(row=1, column=3)
        verbose_spinbox.bind('<Return>', lambda event: validate_spinbox_verbose_input(verbose_spinbox.get()))
        ttk.Label(misc_frame, text="Flip Images:").grid(row=2, column=2, sticky="nsew")
        menu_button = ttk.Menubutton(misc_frame, textvariable=self.flip_var, width=5)
        menu_button.grid(row=2, column=3)
        menu_button.menu = tk.Menu(menu_button, tearoff=0)
        menu_button["menu"] = menu_button.menu
        for item in ["None", "X", "Y", "Z", "XY", "XZ", "YZ", "XYZ"]:
            menu_button.menu.add_radiobutton(label=item, variable=self.flip_var, value=item)

        # Entry for the conf path
        conf_frame = ttk.Frame(misc_frame)
        conf_frame.grid(row=3, column=0, columnspan=4, pady=3, sticky="nsew")
        ttk.Label(conf_frame, text="Conf Path:").grid(row=0, column=0, sticky="nsew", padx=5)
        ttk.Entry(conf_frame, textvariable=self.configuration_path_var,width=25,validate="focusout",validatecommand=lambda: self.validate_folder_path(self.configuration_path_var)).grid(row=0, column=1, pady=3)
        ttk.Button(conf_frame, image=self.folder_icon, command=lambda: self.browse_path(self.configuration_path_var)).grid(row=0, column=2, sticky="e", pady=3, padx=(0, 5))

    def create_algorithm_widgets(self, options_frame):
        algorithm_frame = ttk.LabelFrame(options_frame, text="Algorithm")
        algorithm_frame.grid(row=1, column=1, padx=10, pady=5, sticky="nsew")

        def validate_spinbox_penalty_strength(P):
            # P is the new value of the Spinbox after the edit
            if P == "":
                self.penalty_strength_var.set(0)
                return True
            elif P.isdigit() or self.is_valid_float(P):
                if float(P) > 100:
                    self.penalty_strength_var.set(100)
                    return True
                elif float(P) < 0:
                    self.penalty_strength_var.set(0)
                    return True
                else:
                    return True
            else:
                self.penalty_strength_var.set(0)
                return True
        vcmd_penalty_strength = (self.register(validate_spinbox_penalty_strength), '%P')

        ttk.Label(algorithm_frame, text="Iterations and Subsets:").grid(row=0, column=0, sticky="nsew", padx=(5, 2))
        ttk.Entry(algorithm_frame, textvariable=self.iterations_var,width=29).grid(row=0, column=1, sticky="w")
        ttk.Label(algorithm_frame, text="Format: iterations1:subsets1,iterations2:subsets2,etc...").grid(row=1, column=0, columnspan=2, sticky="nsew", padx=(5, 2))

        opt_proj_pnl_frame = ttk.Frame(algorithm_frame)
        opt_proj_pnl_frame.grid(row=2, column=0, columnspan=2, pady=3, sticky="nsew")
        # Optimizer
        ttk.Label(opt_proj_pnl_frame, text="Optimizer:").grid(row=0, column=0, sticky="nsew", padx=(5, 2))
        opti_menu_button = ttk.Menubutton(opt_proj_pnl_frame, textvariable=self.optimizer_var, width=12)
        opti_menu_button.grid(row=0, column=1, pady=3)
        opti_menu_button.menu = tk.Menu(opti_menu_button, tearoff=0)
        opti_menu_button["menu"] = opti_menu_button.menu
        for item in ["MLEM", "OSL", "DEPIERRO95"]:
            opti_menu_button.menu.add_radiobutton(label=item, variable=self.optimizer_var, value=item)

        # Penalty
        ttk.Label(opt_proj_pnl_frame, text="Penalty:").grid(row=0, column=3, sticky="nsew", padx=(5, 2))
        penalty_menu = ttk.Menubutton(opt_proj_pnl_frame, textvariable=self.penalty_var)
        penalty_menu.grid(row=0, column=4, pady=3, padx=(0, 5))
        penalty_menu.menu = tk.Menu(penalty_menu, tearoff=0)
        penalty_menu["menu"] = penalty_menu.menu
        for item in ["MRF", "MRP"]:
            penalty_menu.menu.add_radiobutton(label=item, variable=self.penalty_var, value=item)

        # Projector
        ttk.Label(opt_proj_pnl_frame, text="Projector:").grid(row=1, column=0, sticky="nsew", padx=(5, 2))
        proj_menu_button = ttk.Menubutton(opt_proj_pnl_frame, textvariable=self.projector_var, width=12)
        proj_menu_button.grid(row=1, column=1)
        proj_menu_button.menu = tk.Menu(proj_menu_button, tearoff=0)
        proj_menu_button["menu"] = proj_menu_button.menu
        for item in ["joseph", "classicSiddon", "incrementalSiddon", "distanceDriven"]:
            proj_menu_button.menu.add_radiobutton(label=item, variable=self.projector_var, value=item)

        # Penalty strength
        ttk.Label(opt_proj_pnl_frame, text="Strength:").grid(row=1, column=3, sticky="nsew", padx=(5, 2))
        spinbox_penalty_strength = ttk.Spinbox(opt_proj_pnl_frame, textvariable=self.penalty_strength_var,from_=0, to=100, increment=0.1, width=6, validate='focusout', validatecommand=vcmd_penalty_strength)
        spinbox_penalty_strength.grid(row=1, column=4, padx=(0, 5))
        spinbox_penalty_strength.bind('<Return>', lambda event: validate_spinbox_penalty_strength(spinbox_penalty_strength.get()))

        # Function to toggle penalty options based on optimizer selection
        def update_penalty_menu():
            # replace with your actual optimizer variable
            optimizer = self.optimizer_var.get()
            if optimizer == "MLEM":
                penalty_menu.config(state="disabled")
                spinbox_penalty_strength.config(state="disabled")
            else:
                penalty_menu.config(state="normal")
                spinbox_penalty_strength.config(state="normal")
                if optimizer == "OSL":
                    penalty_menu.menu.delete(0, 'end')
                    for item in ["MRF", "MRP"]:
                        penalty_menu.menu.add_radiobutton(label=item, variable=self.penalty_var, value=item)
                elif optimizer == "DEPIERRO95":
                    penalty_menu.menu.delete(0, 'end')
                    penalty_menu.menu.add_radiobutton(label="MRF", variable=self.penalty_var, value="MRF")
                    self.penalty_var.set("MRF")
        self.optimizer_var.trace_add("write", lambda *args: update_penalty_menu())
        # Initially call the function to set the initial state of penalty options
        update_penalty_menu()

    def update_convolution_entries(self):
        # Get the current number of convolutions
        num_conv = self.convolution_num_var.get()

        # If the number of convolutions has decreased, remove the extra entries and convolutions
        if num_conv < self.previous_num_conv:
            del self.convolution_type_vars[num_conv:]
            del self.convolution_value_vars[num_conv:]
            del self.convolution_x_var[num_conv:]
            del self.convolution_y_var[num_conv:]
            del self.convolution_sigma_var[num_conv:]
            for i in range(num_conv, self.previous_num_conv):
                for widget in self.convolution_frame.grid_slaves():
                    row = int(widget.grid_info()["row"])
                    col = int(widget.grid_info()["column"])
                    if (row == (i//2)*2 or row == 1+(i//2)*2) and col in range(3, 11):widget.grid_forget()

        # If the number of convolutions has increased, add the extra entries and convolutions
        elif num_conv > self.previous_num_conv:
            for i in range(self.previous_num_conv, num_conv):
                type_var = tk.StringVar()
                type_var.set("psf")
                value_x_var = tk.DoubleVar()
                value_x_var.set(1.0)
                value_y_var = tk.DoubleVar()
                value_y_var.set(1.0)
                value_sigma_var = tk.DoubleVar()
                value_sigma_var.set(3.0)
                value_var = tk.StringVar()
                value_var.set(f"gaussian,{value_x_var.get()},{value_x_var.get()},{value_sigma_var.get()}::{type_var.get()}")
                self.convolution_type_vars.append(type_var)
                self.convolution_value_vars.append(value_var)
                self.convolution_x_var.append(value_x_var)
                self.convolution_y_var.append(value_y_var)
                self.convolution_sigma_var.append(value_sigma_var)

        # Update the entries based on the number of convolutions
        if num_conv == 0:
            ttk.Label(self.convolution_frame, text="Conv 0. Type:").grid(row=0, column=3, padx=(5, 2))
            menu_button = ttk.Menubutton(self.convolution_frame, textvariable=tk.StringVar(), state="disabled", width=10)
            menu_button.grid(row=0, column=4, columnspan=3)
            menu_button.menu = tk.Menu(menu_button, tearoff=0)
            menu_button["menu"] = menu_button.menu
            for item in ["psf", "post", "sieve", "intra", "backward", "forward"]:
                menu_button.menu.add_radiobutton(label=item, variable=self.convolution_type_vars, value=item)
            ttk.Label(self.convolution_frame, text="Size (X,Y,Sig):").grid(row=1, column=3, padx=(5, 2))
            ttk.Spinbox(self.convolution_frame, textvariable=tk.StringVar(), state="disabled", width=3).grid(row=1, column=4)
            ttk.Spinbox(self.convolution_frame, textvariable=tk.StringVar(), state="disabled", width=3).grid(row=1, column=5)
            ttk.Spinbox(self.convolution_frame, textvariable=tk.StringVar(), state="disabled", width=3).grid(row=1, column=6)
        else:
            for i in range(0, num_conv):
                j = 0 if i % 2 == 0 else 1
                ttk.Label(self.convolution_frame, text=f"Conv {i+1} Type:").grid(row=(i//2)*2, column=3+j*4, sticky="nsew", padx=(5, 2))
                menu_button = ttk.Menubutton(self.convolution_frame, textvariable=self.convolution_type_vars[i], width=10)
                menu_button.grid(row=(i//2)*2, column=4+j*4, columnspan=3)
                menu_button.menu = tk.Menu(menu_button, tearoff=0)
                menu_button["menu"] = menu_button.menu
                for item in ["psf", "post", "sieve", "intra", "backward", "forward"]:
                    menu_button.menu.add_radiobutton(label=item, variable=self.convolution_type_vars[i], value=item)
                ttk.Label(self.convolution_frame, text="Size (X,Y,Sig):").grid(row=1+(i//2)*2, column=3+j*4, sticky="nsew", padx=(5, 2))
                ttk.Spinbox(self.convolution_frame, from_=0, to=100, increment=0.1, textvariable=self.convolution_x_var[i], width=3, command=self.update_entries).grid(row=1+(i//2)*2, column=4+j*4, pady=(0, 5))
                ttk.Spinbox(self.convolution_frame, from_=0, to=100, increment=0.1, textvariable=self.convolution_y_var[i], width=3, command=self.update_entries).grid(row=1+(i//2)*2, column=5+j*4, pady=(0, 5))
                ttk.Spinbox(self.convolution_frame, from_=0, to=100, increment=0.1, textvariable=self.convolution_sigma_var[i], width=3, command=self.update_entries).grid(row=1+(i//2)*2, column=6+j*4, pady=(0, 5))

        self.previous_num_conv = num_conv
        self.update_entries()

    def create_convolution_widgets(self, options_frame):
        self.convolution_frame = ttk.LabelFrame(options_frame, text="Convolutions")
        self.convolution_frame.grid(row=2, column=0, columnspan=2,padx=10, pady=5, sticky="nsew")

        def toggle_spinbox():
            if self.convolution_need_bool_var.get() == 0:
                self.conv_spinbox.set(0)
                self.update_convolution_entries()
                self.conv_spinbox.config(state='disabled')
            else:
                self.conv_spinbox.set(1)
                self.update_convolution_entries()
                self.conv_spinbox.config(state='normal')

        def validate_spinbox_input(P):
            # P is the new value of the Spinbox after the edit
            if P == "":
                self.conv_spinbox.set(0)
            elif P.isdigit():
                if int(P) > 10:
                    self.conv_spinbox.set(10)
                elif int(P) < 0:
                    self.conv_spinbox.set(0)
            else:
                self.conv_spinbox.set(0)
                
            self.update_convolution_entries()
            return True
        vcmd = (self.register(validate_spinbox_input), '%P')


        # CheckButton for Convolution
        ttk.Checkbutton(self.convolution_frame, text="Apply Convolutions", variable=self.convolution_need_bool_var,command=toggle_spinbox).grid(row=0, column=0, columnspan=2, sticky="nsew", padx=5)
        # Spinbox for number of convolutions
        ttk.Label(self.convolution_frame, text="Nb of Conv:").grid(row=1, column=0, sticky="nsew", padx=(5, 0), pady=3)
        self.conv_spinbox = ttk.Spinbox(self.convolution_frame, textvariable=self.convolution_num_var, from_=0, to=10,increment=1, width=3, validate='focusout', validatecommand=vcmd, command=lambda: self.update_convolution_entries())
        self.conv_spinbox.grid(row=1, column=1, pady=(0, 5))
        self.conv_spinbox.bind('<Return>', lambda event: validate_spinbox_input(self.conv_spinbox.get()))
        # Separator between these two and convolutions entries
        separator = ttk.Separator(self.convolution_frame, orient="vertical")
        separator.grid(row=0, column=2, rowspan=2, sticky="ns", padx=5, pady=5)
        # Call the function initially to create the initial entries
        self.update_convolution_entries()

    def convert_absolute_to_relative_path(self, abs_path):
        # Assuming your script's directory is the base for relative paths
        rel_path = os.path.relpath(abs_path, self.script_dir)
        return os.path.normpath(rel_path)

    def browse_path(self, path_var):
        initial_dir = path_var.get() if path_var.get() else self.script_dir
        initial_dir = os.path.dirname(os.path.abspath(initial_dir))
        abs_path = filedialog.askopenfilename(initialdir=initial_dir)
        if abs_path:  # Check if abs_path is not empty
            rel_path = self.convert_absolute_to_relative_path(abs_path)
            path_var.set(rel_path)

    def validate_program_path(self, path_var):
        path = path_var.get()
        if path_var == self.sensitivity_path_var and path == "":
            return True
        elif not os.path.isfile(path):
            messagebox.showerror("Invalid File Path", "The specified file path is invalid.")
            return False
        elif path_var == self.main_program_path_var and not path.endswith(".exe"):
            messagebox.showerror("Invalid Main Program Path", "The specified path is not an executable.")
            return False
        elif path_var == self.datafile_path_var and not path.endswith(".Cdh"):
            messagebox.showerror("Invalid DataFile Path", "The specified path is not a Cdh file.")
            return False
        elif path_var == self.sensitivity_path_var and not path.endswith(".hdr"):
            messagebox.showerror("Invalid Sensitivity Path", "The specified file path is not a hdr file.")
            return False
        return True
    
    def validate_folder_path(self, path_var):
        path = path_var.get()
        if path_var == self.configuration_path_var and path == "":
            return True
        elif not os.path.isdir(path):
            messagebox.showerror("Invalid Folder Path", "The specified folder path is invalid.")
            return False
        return True
    
    def is_valid_float(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def generate_script(self, save = True, info = True):
        # Retrieve user inputs
        mpi_enabled = self.mpi_bool_var.get()
        mpi_threads = self.mpi_threads_var.get()
        verbose_level = self.verbose_level_var.get()
        last_iter = self.last_iter_bool_var.get()
        flip = self.flip_var.get()
        stats_need = self.stats_need_bool_var.get()

        main_program_path = self.main_program_path_var.get()
        datafile = self.datafile_path_var.get()
        output_path = self.output_path_var.get()
        sensitivity_path = self.sensitivity_path_var.get()
        configuration_path = self.configuration_path_var.get()

        voxel_number = self.voxel_number_var.get()
        voxel_size = self.voxel_size_var.get()
        fov_size = self.fov_size_var.get()
        offset = self.offset_var.get()
        voxel_size_x = self.voxel_size_x_var.get()
        voxel_size_y = self.voxel_size_y_var.get()
        voxel_size_z = self.voxel_size_z_var.get()
        voxel_number_x = self.voxel_number_x_var.get()
        voxel_number_y = self.voxel_number_y_var.get()
        voxel_number_z = self.voxel_number_z_var.get()
        fov_size_x = self.fov_size_x_var.get()
        fov_size_y = self.fov_size_y_var.get()
        fov_size_z = self.fov_size_z_var.get()
        offset_x = self.offset_x_var.get()
        offset_y = self.offset_y_var.get()
        offset_z = self.offset_z_var.get()

        iterations = self.iterations_var.get()
        optimizer = self.optimizer_var.get()
        projector = self.projector_var.get()
        penalty = self.penalty_var.get()
        penalty_strength = self.penalty_strength_var.get()
        convolution_need = self.convolution_need_bool_var.get()
        convolution_num = self.convolution_num_var.get()
        convolution_types = [var.get() for var in self.convolution_type_vars]
        convolution_values = [var for var in self.convolution_value_vars]
        convolution_x = [var.get() for var in self.convolution_x_var]
        convolution_y = [var.get() for var in self.convolution_y_var]
        convolution_sigma = [var.get() for var in self.convolution_sigma_var]

       
        # Start writing the batch script
        script_content = "@echo off\n"
        script_content += ":: Batch script generated by Python GUI to run the CASToR program\n"
        script_content += ":: To get help run the main programm in the comman-line withou any argument or with 'h', '-help' or '--help' options\n\n"

        # Set Command Line Options
        script_content += ":::::::::::::::::::::::::::::\n"
        script_content += ":: Set Command Line Options\n"
        script_content += ":::::::::::::::::::::::::::::\n\n"

        # Set MPI stuff and recon program
        if mpi_enabled:
            script_content += f"@REM MPI enabled with {mpi_threads} threads\n"
            if mpi_threads > 0:
                script_content += "set mpi_exe=mpiexec.exe -n %d\n" % mpi_threads
            if mpi_threads == 0:
                script_content += "set mpi_exe=mpiexec.exe\n"
        else:
            script_content += "@REM MPI disabled\n"
            script_content += "set mpi_exe=\n"

        # Set the main program path
        script_content += f"set recon_exe={main_program_path}\n\n"

        # Set verbose level, threads, stats need, last iter, flip
        script_content += f"set verbose= -vb {verbose_level}\n"
        script_content += f"set threads= -th {mpi_threads}\n\n"
        if stats_need:
            script_content += "set stats= -opti-stat\n"
        else: 
            script_content += "set stats=\n"
        if last_iter:
            script_content += "set last_it= -oit -1\n"
        else:
            script_content += "set last_it=\n"
        if flip != "None":
            script_content += f"set flip_out= -flip-out {flip}\n"
        else:
            script_content += "set flip_out=\n\n"
        
        # Set the datafile path
        script_content += f"set datafile= -df {datafile}\n"
        # Set the output path
        script_content += f"set output= -dout {output_path}\n"
        # Set the sensitivity path
        if sensitivity_path != "":
            script_content += f"set sensitivity= -sens {sensitivity_path}\n"
        else:
            script_content += "set sensitivity=\n"
        # Set the configuration path
        if configuration_path != "":
            script_content += f"set configuration= -conf {configuration_path}\n\n"
        else:
            script_content += "set configuration=\n\n"

        # Set the voxel number, size, fov size, offset
        script_content += f"set voxel_number= -dim {voxel_number_x},{voxel_number_y},{voxel_number_z}\n"
        script_content += f"set voxel_size= -vox {voxel_size_x},{voxel_size_y},{voxel_size_z}\n"
        script_content += f"set offset= -off {offset_x},{offset_y},{offset_z}\n\n"

        # Set the iterations, optimizer, projector, penalty, penalty strength
        script_content += f"set iterations= -it {iterations}\n"
        script_content += f"set optimizer= -opti {optimizer}\n"
        script_content += f"set projector= -proj {projector}\n"
        if optimizer == "MLEM":
            script_content += "set penalty=\n"
            script_content += "set penalty_strength=\n"
        else:
            script_content += f"set penalty= -pnlt {penalty}\n"
            script_content += f"set penalty_strength= -pnlt-beta {penalty_strength}\n\n"
        
        # Set the convolution options
        if convolution_need:
            script_content += f"@REM Set {convolution_num} convolution(s)\n"
            for i in range(convolution_num):
                script_content += f"set psf_{i+1}= -conv {convolution_values[i]}\n\n"

        # Run the main program
        script_content += ":::::::::::::::::::::::::::::\n"
        script_content += ":: Launch the reconstruction\n"
        script_content += ":::::::::::::::::::::::::::::\n\n"

        script_content += "echo ==========================================================\n"
        script_content += "echo Reconstruction is going on. Should take several minutes\n"
        script_content += "echo ==========================================================\n\n"

        script_content += f"%mpi_exe% %recon_exe% %verbose% %threads% %stats% %last_it% %flip_out% %datafile% %output% %sensitivity% %configuration% %voxel_number% %voxel_size% %offset% %iterations% %optimizer% %projector% %penalty% %penalty_strength% "
        if convolution_need:
            for i in range(convolution_num):
                script_content += f"%psf_{i+1}% "
        script_content += "\n\n"

        script_content += "echo ==========================================================\n"
        script_content += "echo Reconstruction is finished!\n"
        script_content += "echo ==========================================================\n\n"
        # Define the path for the script
        self.script_name = "run_castor_python_script.bat"
        self.script_path = os.path.join(self.script_dir, self.script_name)
        
        if save:
            # Write and save the script to the file
            with open(self.script_path, "w") as file:
                file.write(script_content)
        
        # Display the generated script in a message box
        if info:
            with open(self.script_path, "r") as file:
                script = file.read()
            messagebox.showinfo("Generated Batch Script", script)

    def run_castor_program(self):
        # Run CASToR with the generated script
        self.generate_script(save=False, info=True)
        print(f'Running the CASToR program with the generated script: \"{self.script_path}\"')
        command = f"powershell -NoExit -Command \"& '{self.script_path}'\""
        os.system(f"start {command}")

    def print_test_all_variables(self):
        self.update_entries()

        print("All Variables:")
        print("MPI:", self.mpi_bool_var.get())
        print("MPI Threads:", self.mpi_threads_var.get())
        print("Verbose Level:", self.verbose_level_var.get())
        print("Last Iteration:", self.last_iter_bool_var.get())
        print("Flip:", self.flip_var.get())
        print("Stats Need:", self.stats_need_bool_var.get())

        print("Main Program Path:", self.main_program_path_var.get())
        print("Datafile Path:", self.datafile_path_var.get())
        print("Output Path:", self.output_path_var.get())
        print("Sensitivity Path:", self.sensitivity_path_var.get())
        print("Configuration Path:", self.configuration_path_var.get())

        print("Voxel Number:", self.voxel_number_var.get())
        print("Voxel Size:", self.voxel_size_var.get())
        print("FOV Size:", self.fov_size_var.get())
        print("Offset:", self.offset_var.get())
        print("Voxel Size X:", self.voxel_size_x_var.get())
        print("Voxel Size Y:", self.voxel_size_y_var.get())
        print("Voxel Size Z:", self.voxel_size_z_var.get())
        print("Voxel Number X:", self.voxel_number_x_var.get())
        print("Voxel Number Y:", self.voxel_number_y_var.get())
        print("Voxel Number Z:", self.voxel_number_z_var.get())
        print("FOV Size X:", self.fov_size_x_var.get())
        print("FOV Size Y:", self.fov_size_y_var.get())
        print("FOV Size Z:", self.fov_size_z_var.get())
        print("Offset X:", self.offset_x_var.get())
        print("Offset Y:", self.offset_y_var.get())
        print("Offset Z:", self.offset_z_var.get())

        print("Iterations:", self.iterations_var.get())
        print("Optimizer:", self.optimizer_var.get())
        print("Projector:", self.projector_var.get())
        print("Penalty:", self.penalty_var.get())
        print("Penalty Strength:", self.penalty_strength_var.get())
        print("Convolution Need:", self.convolution_need_bool_var.get())
        print("Convolution Num:", self.convolution_num_var.get())
        print("Convolution Types:", [var.get()for var in self.convolution_type_vars])
        print("Convolution Values:", [var for var in self.convolution_value_vars])
        print("Convolution X:", [var.get() for var in self.convolution_x_var])
        print("Convolution Y:", [var.get() for var in self.convolution_y_var])
        print("Convolution Sigma:", [var.get()for var in self.convolution_sigma_var])

    def show_about(self):
        help_window = tk.Toplevel()
        help_window.title("About Me")
        help_window.geometry("1050x350")
        help_window.resizable(False, False)

        help_text = """
        This is a GUI application to generate a batch script for the CASToR program.
        The CASToR program is a reconstruction software for Positron Emission Tomography (PET) and Computed Tomography (CT) images.
        The batch script generated by this application can be used to creates scripts and run the CASToR program with the specified parameters.

        The application has several input fields and options to configure the parameters for the CASToR program:
        The user can set the main program path, datafile path, output path, sensitivity path, configuration path, and other parameters.
        The user can also set the voxel number, size, FOV size, offset, iterations, optimizer, projector, penalty, and convolution options.

        The application provides a 'Generate Script' button to generate the batch script based on the user inputs.
        The user can save the generated script to a file and run the CASToR program with the generated script.

        The application also provides a 'Run CASToR Program' button to directly run the CASToR program with the generated script.

        The user can also print all the variables set in the application to check the values before generating the script.

        For more information about the CASToR program and its parameters, please refer to the CASToR documentation (https://castor-project.org/).
        """

        help_label = tk.Label(help_window, text=help_text, justify="left", font=("Arial", 12))
        help_label.pack(padx=5, pady=10, fill="both", expand=True)

    def open_from_file(self):
        file_path = filedialog.askopenfilename(initialdir=self.script_dir)
        if file_path:
            with open(file_path, "r") as file:
                lines = file.readlines()
                del self.convolution_value_vars[:], self.convolution_type_vars[:], self.convolution_x_var[:], self.convolution_y_var[:], self.convolution_sigma_var[:]
                self.convolution_num_var.set(0)
                for line in lines:
                    print(line)
                    if "set mpi_exe=" in line:
                        mpi_exe = line.split("set mpi_exe=")[1]
                        if mpi_exe == "\n":
                            self.mpi_bool_var.set(False)
                        else:
                            self.mpi_bool_var.set(True)
                    if "set threads=" in line:
                        mpi_threads = line.split("set threads=")[1].split()[1]
                        self.mpi_threads_var.set(mpi_threads)
                    if "set verbose=" in line:
                        verbose = line.split("set verbose=")[1].split()[1]
                        self.verbose_level_var.set(verbose)
                    if "set last_it=" in line:
                        last_it = line.split("set last_it=")[1]
                        if last_it == "\n":
                            self.last_iter_bool_var.set(False)
                        else:
                            self.last_iter_bool_var.set(True)
                    if "set flip_out=" in line:
                        flip_out = line.split("set flip_out=")[1]
                        if flip_out == "\n":
                            self.flip_var.set("None")
                        else:
                            self.flip_var.set(flip_out.split()[1])
                    if "set stats=" in line:
                        stats = line.split("set stats=")[1]
                        if stats == "\n":
                            self.stats_need_bool_var.set(False)
                        else:
                            self.stats_need_bool_var.set(True)
                    if "set recon_exe=" in line:
                        recon_exe = line.split("set recon_exe=")[1]
                        self.main_program_path_var.set(recon_exe)
                    if "set datafile=" in line:
                        datafile = line.split("set datafile=")[1].split()[1]
                        self.datafile_path_var.set(datafile)
                    if "set output=" in line:
                        output = line.split("set output=")[1].split()[1]
                        self.output_path_var.set(output)
                    if "set sensitivity=" in line:
                        sensitivity = line.split("set sensitivity=")[1]
                        if sensitivity == "\n":
                            self.sensitivity_path_var.set("")
                        else:
                            self.sensitivity_path_var.set(sensitivity.split()[1])
                    if "set configuration=" in line:
                        configuration = line.split("set configuration=")[1]
                        if configuration == "\n":
                            self.configuration_path_var.set("")
                        else:
                            self.configuration_path_var.set(configuration.split()[1])
                    if "set voxel_number=" in line:
                        voxel_number = line.split("set voxel_number=")[1].split()[1]
                        self.voxel_number_var.set(voxel_number)
                        # Update the voxel number x, y, z
                        self.voxel_number_x_var.set(voxel_number.split(",")[0])
                        self.voxel_number_y_var.set(voxel_number.split(",")[1])
                        self.voxel_number_z_var.set(voxel_number.split(",")[2])
                    if "set voxel_size=" in line:
                        voxel_size = line.split("set voxel_size=")[1].split()[1]
                        self.voxel_size_var.set(voxel_size)
                        # Update the voxel size x, y, z
                        self.voxel_size_x_var.set(voxel_size.split(",")[0])
                        self.voxel_size_y_var.set(voxel_size.split(",")[1])
                        self.voxel_size_z_var.set(voxel_size.split(",")[2])
                    if "set offset=" in line:
                        offset = line.split("set offset=")[1].split()[1]
                        self.offset_var.set(offset)
                        # Update the offset x, y, z
                        self.offset_x_var.set(offset.split(",")[0])
                        self.offset_y_var.set(offset.split(",")[1])
                        self.offset_z_var.set(offset.split(",")[2])
                    if "set iterations=" in line:
                        iterations = line.split("set iterations=")[1].split()[1]
                        self.iterations_var.set(iterations)
                    if "set optimizer=" in line:
                        optimizer = line.split("set optimizer=")[1].split()[1]
                        self.optimizer_var.set(optimizer)
                    if "set projector=" in line:
                        projector = line.split("set projector=")[1].split()[1]
                        self.projector_var.set(projector)
                    if "set penalty=" in line:
                        penalty = line.split("set penalty=")[1]
                        if penalty == "\n":
                            self.penalty_var.set("MRF")
                        else:
                            self.penalty_var.set(penalty.split()[1])
                    if "set penalty_strength=" in line:
                        penalty_strength = line.split("set penalty_strength=")[1]
                        if penalty_strength == "\n":
                            self.penalty_strength_var.set(0.5)
                        else:
                            self.penalty_strength_var.set(penalty_strength.split()[1])
                    if "set psf_" in line:
                        psf = line.split("=")[1].strip().split()[1]  # Split the line at "=" and remove leading/trailing whitespace
                        print(psf)
                        if psf == "\n":
                            self.convolution_need_bool_var.set(False)
                            self.conv_spinbox.config(state='disabled')
                        else:
                            self.convolution_need_bool_var.set(True)
                            self.conv_spinbox.config(state='normal')
                            self.convolution_num_var.set(self.convolution_num_var.get() + 1)
                            # self.convolution_value_vars.append(psf)
                            psf_values = psf.split("::")[0].split(",")  # Split the psf values at "::" and ","
                            self.convolution_type_vars.append(tk.StringVar(value=psf.split("::")[1]))
                            self.convolution_x_var.append(tk.DoubleVar(value=float(psf_values[1])))
                            self.convolution_y_var.append(tk.DoubleVar(value=float(psf_values[2])))
                            self.convolution_sigma_var.append(tk.DoubleVar(value=float(psf_values[3])))
                        self.update_convolution_entries()
        
        # Get from the file the values of the variables

    
if __name__ == "__main__":
    app = BatchScriptGenerator()
    app.mainloop()
