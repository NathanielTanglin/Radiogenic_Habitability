import os
import numpy as np

class Parameter_Sweep:
    vspace_file_parse = '# dummy file workaround for BigPlanet.\nsDestFolder {directory_name}\nsTrialName {trial_name}\n\n{bodies}\n{parameters}\nsPrimaryFile {primary_file}'

    # The parameter sweep creates a R x C matrix with as many columns
    # (C) as varied parameters and as many rows (R) as permutations over
    # the range of variations.

    # All the data is indexed by column. Each column will have a corresponding
    # input file and variations data.

    def __init__(self, **kwargs):
        if 'trial_name' in kwargs:
            self.trial_name = kwargs['trial_name']
        else:
            self.trial_name = 'run'

        self.input_options = []
        self.ranges = []
        self.names = []
        self.input_file_paths = []
        self.unvaried_input_files = []

        # The keys of the input files dictionary map to the variations data for each varied input option.
        for path in kwargs['paths']:
            # Splices out the file name.
            file_name = path.split(os.sep)[-1]

            # Checks that the file name is a key in the arguments dictionary (i.e. if that specific input file
            # contains a varied parameter).
            if file_name in kwargs:
                # Populates a set of parallel arrays.
                self.names += kwargs[file_name]['names']
                self.ranges += kwargs[file_name]['ranges']

                for input_options in kwargs[file_name]['input_options']:
                    self.input_options.append(input_options)
                    self.input_file_paths.append(path)
            else:
                self.unvaried_input_files.append(path)

    def _num_combinations(self, n, p):
        answer = 1
        
        for i in range(n, len(p)):
            answer *= len(p[i])

        return answer
    
    def _matrix_component(self, row, column, matrix, return_k=False):
        if row+1 > self._num_combinations(0, matrix):
            raise Exception('Row out of bounds of matrix.')

        k = 0

        if column+1 == len(matrix):
            k = row % len(matrix[column])
        else:
            k = int(int(np.ceil((row+1)/self._num_combinations(column+1, matrix))-1) % len(matrix[column]))

        if return_k:
            return (k, matrix[column][k])
        else:
            return matrix[column][k]

    def _inject_input_file(self, input_file_path, injection):
        with open(input_file_path, 'r') as input_file:
            contents = input_file.read()
            (start, sep, end) = contents.partition('saOutputOrder')

            unwanted_strings = ('\n', ' ')

            # Remove all '\n' and ' ' strings from end of the 'start' partition to guarantee neat spacing.
            while start.endswith(unwanted_strings):
                for suffix in unwanted_strings:
                    start = start.removesuffix(suffix)

            while end.startswith(unwanted_strings):
                for prefix in unwanted_strings:
                    end = end.removeprefix(prefix)

            # Remove all '\n' and ' ' strings from end of the 'start' partition to guarantee neat spacing.
            while injection.endswith(unwanted_strings):
                for suffix in unwanted_strings:
                    injection = injection.removesuffix(suffix)

            while injection.startswith(unwanted_strings):
                for prefix in unwanted_strings:
                    injection = injection.removeprefix(prefix)
        
            if sep != '':
                output = start + '\n\n' + injection + '\n\n' + sep + ' ' + end
            else:
                output = start + '\n\n' + injection + sep + ' ' + end

            return output
    
    def generate_input_files(self, directory_path = 'Parameter_Sweep'):
        directory_name = directory_path.split(os.sep)[-1]
        parent_directory = directory_path.removesuffix(directory_name)

        # Creates the directory in which the input files will be placed.
        if directory_name not in os.listdir(parent_directory):
            os.mkdir(directory_path)

        for row in range(self._num_combinations(0, self.ranges)):
            injections = {}
            
            # The file name for the subdirectory for this specific parameter combination.
            parameter_file_identifier = ''

            for column in range(len(self.ranges)):
                (k, abundance) = self._matrix_component(row, column, self.ranges, return_k = True)
                
                parameter_file_identifier += '_{parameter_name}{n}'.format(parameter_name = self.names[column], n = k)

                path = self.input_file_paths[column]
                contents = open(path).read()

                if path not in injections:
                    injections[path] = ''

                injections[path] += '# Varied parameter {n}\n'.format(n=column+1)

                for option in self.input_options[column]:
                    #new_string = '{option} -{abundance}\n'.format(option=option, abundance=abundance)
                    #if contents.find(option) != -1:

                    unit_character = ''

                    if option[0] == '-':
                        option = option.removeprefix('-')
                        unit_character = '-'

                    injections[path] += '{option} {abundance}\n'.format(option=option, abundance=unit_character+str(abundance))
                    #else:
                        #contents.replace(option, '')

                injections[path] += '\n'

            # Subdirectory for the specific run. This is where the modified input files go.
            run_name = '{trial}{parameter_file_identifier}'.format(trial = self.trial_name, parameter_file_identifier=parameter_file_identifier)
            run_directory = os.path.join(directory_path, run_name)

            # Makes the new directory.
            if run_name not in os.listdir(directory_path):
                os.mkdir(run_directory)

            for (file_path, injection) in injections.items():
                file_name = file_path.split(os.sep)[-1]

                body_file_path = os.path.join(run_directory, file_name)

                with open(body_file_path, 'w', encoding = 'utf-8') as input_file:
                    contents = self._inject_input_file(file_path, injection)

                    input_file.write(contents)

            for file_path in self.unvaried_input_files:
                contents = ''

                with open(file_path, 'r') as input_file:
                    contents = input_file.read()

                file_name = file_path.split(os.sep)[-1]
                body_file_path = os.path.join(run_directory, file_name)

                with open(body_file_path, 'w') as input_file:
                    input_file.write(contents)

        vspace_file_path = os.path.join(directory_path, os.pardir, 'vspace.in')

        varied_file_paths = self.input_file_paths
        unvaried_file_paths = self.unvaried_input_files

        bodies = ''
        primary_file = ''

        for input_file_path in set(varied_file_paths + unvaried_file_paths):
            with open(input_file_path) as input_file:
                contents = input_file.read()
                input_file_name = input_file_path.split(os.sep)[-1]

                if contents.find('saBodyFiles') != -1:
                    primary_file = input_file_name
                
                if input_file_path in varied_file_paths:
                    bodies += 'sBodyFile ' + input_file_name + '\n'

        parameters = ''

        for (n, name) in enumerate(self.names):
            parameters += 'PARAMETER_{n} [0, 0, n0] {name}\n'.format(n=n, name=name)

        with open(vspace_file_path, 'w', encoding = 'utf-8') as vpl_file:
            vpl_file.write(self.vspace_file_parse.format(directory_name=directory_name, bodies=bodies, parameters=parameters, trial_name=self.trial_name, primary_file=primary_file))