#!/usr/bin/env python

############### ROBOT ASSEMBLY CLASS ##################################
#!/usr/bin/env python

# For automating part
import synbiopython.lab_automation as lab

''' This part of the design module is used for robot instructions - primarily for Flowbot one

class RobotAssembly
Input:
    :param Pandas_dataframe_PCR
    :param F_primers: list
    :param R_primers: list
    :param:Templates:list

    Robot_output:
    :param  self.source_plateF_primers, self.source_plateR_primers, self.source_plateTemplates, self.source_platePCRmix = sourceplate objets
    :param  self.picklist
-----


HELPER FUNCTIONS: 
-----------------
well_keys_96
MakeVirtualPlatesFromDF
PicklistFromPlates

'''


# Helper functions

def well_keys_96(row = True):
    '''
    If true it generates keys for a 96 well plate by row.
    else it does it by column
    '''
    key_list = 'ABCDEFGH'

    if row == True:
        wells_keys = []
        for i in range(0,len(key_list)):
            for j in range(1,13):
                wells_keys.append(key_list[i]+str(j))

    else:
        wells_keys = []
        for i in range(1,13):
            for j in range(0,len(key_list)):
                wells_keys.append(str(key_list[j])+str(i))

    return wells_keys

well_keys = well_keys_96()

def make_virtual_plates_fromDF(f_primers:list, r_primers:list, templates:list, Dataframe_with_PCR_contents):
    '''
    THis function can make virtual plates from lists of primers and templates.
    The Pandas DataFrame is used to calculate how much pcr and h2o is needed for the reactions.
    '''

    #initializing
    source_plate1 = lab.Plate96(name="1") # Forward primers
    source_plate2 = lab.Plate96(name="2") # Reverse primers
    source_plate3 = lab.Plate96(name="3") # Templates
    source_plate4 = lab.Plate2x4(name="4") # PCR mix and H20
    well_keys = well_keys_96()
    volume_needed = len(Dataframe_with_PCR_contents)+1

    # Adding each primer/template to the virtual plates here a volume of 5 is arbitrarily chosen
    for i in range(0,len(f_primers)):
        source_plate1[well_keys[i]].add_content({f_primers[i]:1},volume= 5.0)
    for j in range(0,len(r_primers)):
        source_plate2[well_keys[j]].add_content({r_primers[j]:1},volume= 5.0)
    for k in range(0,len(templates)):
        source_plate3[well_keys[k]].add_content({templates[k]:1},volume= 5.0 )

    # Adding MM and H2O to the final plate
    source_plate4["A1"].add_content({'PCR_2x_mix':10},volume = volume_needed*10)
    source_plate4["A2"].add_content({'H2O_MQ':10},volume = volume_needed*7)

    return source_plate1, source_plate2, source_plate3, source_plate4


def picklist_from_plates(F_primer_plate, R_primer_plate, Templates_plate, MM_H20_plate, PCR_dataframe):
    '''This function can generate picklist from virtual plates and pandas dataframe with PCR components'''

    FORWARD_PRIMERS = F_primer_plate.to_pandas_dataframe()
    REVERSE_PRIMERS = R_primer_plate.to_pandas_dataframe()
    TEMPLATES = Templates_plate.to_pandas_dataframe()
    MM_H2O = MM_H20_plate.to_pandas_dataframe()

    destination_plate = lab.Plate96(name="5")
    picklist = lab.PickList()
    counter = 0

    # FIRST WE CHECK WICH COMPONENTS WE NEED FOR THE PCRs
    for index,row in PCR_dataframe.iterrows():
        templates = row['Template']
        f_primers = row['forward_primer']
        r_primers = row['reverse_primer']

        # FORWARD PRIMERS
        # we see if we can find the primer we need in the primer plate and if we find it we will add it to the picklist
        for index1, row1 in FORWARD_PRIMERS.iterrows():
            for k, v in  row1['content']['quantities'].items():
                if k == f_primers:
                    transfer_forward_primers = lab.Transfer(F_primer_plate.wells[index1], destination_plate.wells[well_keys[counter]], 1 )
                    picklist.add_transfer(transfer = transfer_forward_primers)

        # REVERSE PRIMERS
        # we see if we can find the primer we need in the primer plate
        for index2, row2 in REVERSE_PRIMERS.iterrows():
            for k, v in row2['content']['quantities'].items():
                if k == r_primers:
                    transfer_reverse_primers = lab.Transfer(R_primer_plate.wells[index2], destination_plate.wells[well_keys[counter]], 1 )
                    picklist.add_transfer(transfer = transfer_reverse_primers)

        # TEMPLATE
        # we see if we can find the primer we need in the plate
        for index3, row3 in TEMPLATES.iterrows():
            for k, v in row3['content']['quantities'].items():
                if k == templates:
                    transfer_template = lab.Transfer(Templates_plate.wells[index3], destination_plate.wells[well_keys[counter]], 1 )
                    picklist.add_transfer(transfer = transfer_template)

        # Master mix
        transfer_MM = lab.Transfer(MM_H20_plate.wells['A1'], destination_plate.wells[well_keys[counter]], 10 )
        picklist.add_transfer(transfer = transfer_MM)
        #H20
        transfer_H20 = lab.Transfer(MM_H20_plate.wells['A2'], destination_plate.wells[well_keys[counter]], 7 )
        picklist.add_transfer(transfer = transfer_H20)

        # TO MAKE SURE WE ADD TO THE CORRECT WELL
        counter += 1

    return picklist


class RobotAssembly:
    """Class to generate instructions for robot on demand. It can be initialized in two ways:
    1. Use a design assembly object
    2. Pandas dataframe with PCRs and Pandas dataframe with Forward and Reverse primers, and Template

    Input:
    :param DesignAssembly_object
    :param Pandas_dataframe_PCR
    :param F_primers: list
    :param R_primers: list
    :param:Templates:list

    Robot_output:
    :param  self.source_plateF_primers, self.source_plateR_primers, self.source_plateTemplates, self.source_platePCRmix = sourceplate objets
    :param  self.picklist
    """

    def __init__(self, Pandas_dataframe_PCR, F_primers:list, R_primers:list, Templates:list):

        ###  1.INITIALIZING ##
        ######################## Assembly from plates ############################
        if Pandas_dataframe_PCR is not None and F_primers is not None and R_primers is not None and Templates is not None:
            # Initializing instances
            self.pandas_PCR = Pandas_dataframe_PCR
            self.forward_primers = F_primers
            self.reverse_primers = R_primers
            self.templates = Templates

            #Virtual plates
            self.source_plateF_primers, self.source_plateR_primers, self.source_plateTemplates, self.source_platePCRmix = make_virtual_plates_fromDF( self.forward_primers, self.reverse_primers, self.templates, self.pandas_PCR)

            # Generating Picklist from virtual plates
            self.picklist = picklist_from_plates(self.source_plateF_primers, self.source_plateR_primers, self.source_plateTemplates, self.source_platePCRmix, self.pandas_PCR)

        else:
            print('Remember to put in all the neccesarry components')
    ##########################################################
    ###################  CLASS METHODS  ######################

    def PlatesToExcelFile(self):
        '''Returns an excel file of the plate setup that needs to be made before the flowbot can operate'''
        df1 = self.source_plateF_primers.to_pandas_dataframe()
        df2 = self.source_plateR_primers.to_pandas_dataframe()
        df3 = self.source_plateTemplates.to_pandas_dataframe()
        df4 = self.source_platePCRmix.to_pandas_dataframe()

        with pd.ExcelWriter('Plate_instructions.xlsx') as writer:
            return df1.to_excel(writer, sheet_name='Forward_primer_wells'), df2.to_excel(writer, sheet_name='Reverse_primer_wells'), df3.to_excel(writer, sheet_name='Template_wells'),df4.to_excel(writer, sheet_name='PCR_MIX_H20')


    def print_well_df_to_string(self):
        ''' Prints Pandas dataframe in string format'''

        df1 = self.source_plateF_primers.to_pandas_dataframe()
        df2 = self.source_plateR_primers.to_pandas_dataframe()
        df3 = self.source_plateTemplates.to_pandas_dataframe()
        df4 = self.source_platePCRmix.to_pandas_dataframe()

        return print(' ###Forward primers: ', df1.to_csv(sep=' ', index=False, header=True),'\n', '###Reverse primers: ', df2.to_csv(sep=' ', index=False, header=True),'\n','###Templates: ',df3.to_csv(sep=' ', index=False, header=True),  df4.to_csv(sep=' ', index=False, header=True))

    def FlowbotInstructionsToCSV(self):
        '''Function that prints flowbot instructions to csv format'''

        with open('Flowbot_instructions.csv','w') as outfile:

            return print('source, target, volume','\n',*self.picklist.to_flowbot_instructions(), sep='',end='\n', file = outfile)
