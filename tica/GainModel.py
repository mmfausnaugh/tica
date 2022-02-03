#
#Copyright (c) 2021 Michael Fausnaugh
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the Software), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.



class GainModel(object):
    """
    class properties--hard code the gain per cam/ccd/slice
    """
    gains ={'cam1':{'ccd1':{"A":5.22,
                            "B":5.21,
                            "C":5.21,
                            "D":5.26},
                    'ccd2':{"A":5.27,
                            "B":5.14,
                            "C":5.11,
                            "D":5.19},
                    'ccd3':{"A":5.32,
                            "B":5.23,
                            "C":5.2 ,
                            "D":5.24},
                    'ccd4':{"A":5.35,
                            "B":5.22, 
                            "C":5.22, 
                            "D":5.18} 
                },
            'cam2':{'ccd1':{"A":5.31,
                            "B":5.24,
                            "C":5.23,
                            "D":5.24},
                    'ccd2':{"A":5.22,
                            "B":5.28,
                            "C":5.32,
                            "D":5.22},
                    'ccd3':{"A":5.28,
                            "B":5.26,
                            "C":5.25,
                            "D":5.21},
                    'ccd4':{"A":5.33,
                            "B":5.2 ,
                            "C":5.3 ,
                            "D":5.23}
                } ,
            'cam3':{'ccd1':{"A":5.25,
                            "B":5.24,
                            "C":5.24,
                            "D":5.24},
                    'ccd2':{"A":5.36,
                            "B":5.26,
                            "C":5.36,
                            "D":5.29},
                    'ccd3':{"A":5.26,
                            "B":5.16,
                            "C":5.22,
                            "D":5.24},
                    'ccd4':{"A":5.17,
                            "B":5.15,
                            "C":5.15,
                            "D":5.17}
                } ,
            'cam4':{'ccd1':{"A":5.26,
                            "B":5.18,
                            "C":5.2,
                            "D":5.2},
                    'ccd2':{"A":5.25,
                            "B":5.17,
                            "C":5.12,
                            "D":5.2},
                    'ccd3':{"A":5.3,
                            "B":5.18,
                            "C":5.27,
                            "D":5.19,},
                    'ccd4':{"A":5.24,
                            "B":5.12,
                            "C":5.16,
                            "D":5.16}
                } 
        }
    def __init__(self,modelfile):
        pass
        
