//
//  idxAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/18/17.
//
//


#include "idxAnalysis.hpp"

using namespace llvm;

//#define IDX_DEBUG

namespace idxAnalysis {
    
    char IndexAnalysis::ID = 0;
    static RegisterPass<IndexAnalysis> X("idxAnalysis", "Array index analysis Pass", false, false);
    
    IndexAnalysis::IndexAnalysis() : FunctionPass(ID), arrayName(), arrayExpression() {}
    
    /* Pushing all instructions that related to the index expression of a load/store instruction into a stack (vector).
       This function will be called by function computeExpression() */
    std::vector<Instruction*> IndexAnalysis::instStackInit(Instruction* inst) {
        std::vector<Instruction*> instStack;
        instStack.push_back(inst);
        
        /* push all the instructions related to the target instruction into stack */
        unsigned long start = 0;
        while(instStack.size() != start) {
            
            if (instStack[start]->getNumOperands() == 1) {
                if (isa<Instruction>(instStack[start]->getOperand(0))) {
                    Instruction* pinst = dyn_cast<Instruction>(instStack[start]->getOperand(0));
                    if (std::find(instStack.begin(), instStack.end(), pinst) == instStack.end())     {
                        if (!isa<PHINode>(pinst)) {
                            instStack.push_back(pinst);
                        }
                    }
                }
            }
            
            if (instStack[start]->getNumOperands() == 2) {
                if (isa<Instruction>(instStack[start]->getOperand(0))) {
                    Instruction* pinst = dyn_cast<Instruction>(instStack[start]->getOperand(0));
                    if (std::find(instStack.begin(), instStack.end(), pinst) == instStack.end()) {
                        if (!isa<PHINode>(pinst)) {
                            instStack.push_back(pinst);
                        }
                    }
                }
                if (isa<Instruction>(instStack[start]->getOperand(1))) {
                    Instruction* pinst = dyn_cast<Instruction>(instStack[start]->getOperand(1));
                    if (std::find(instStack.begin(), instStack.end(), pinst) == instStack.end()) {
                        if (!isa<PHINode>(pinst)) {
                            instStack.push_back(pinst);
                        }
                    }
                }
            }
            
            start++;
            if (start > 100) {
#ifdef IDX_DEBUG
                errs() << "ERROR: Expression should not be that long\n";
#endif
                break;
            }
        }
        return instStack;
    }
    
    
    /* Extracting array index expression from load/store instruction 
       1. Pushing all related instructions into stack by calling instStackInit()
       2. Scanning instruction stack to get a prefix expression 
       3. Convert prefix expression to infix expression */
    std::string IndexAnalysis::computeExpression(Instruction *inst) {
        
        std::vector<Instruction*> instStack = instStackInit(inst);
        
        if (instStack.size() == 0) {
            return "Warning: No expression";
        }
        
        /* iterate over the instructions stacked to get the prefix expression */
        std::vector<std::string> expr;
        
        for (std::vector<Instruction*>::iterator it = instStack.begin(), eit = instStack.end(); it != eit; ++it) {
            
            Instruction* pInst = *it;
            
#ifdef IDX_DEBUG
            pInst->dump();
            errs() << expr.size() << "\n";
#endif
            
            if (pInst->isCast()) {
                //stack[i]->getOperand(0)->dump();
                if (isa<PHINode>(pInst->getOperand(0))) {
                    if (pInst->getOperand(0)->hasName()) {
                        expr.push_back(pInst->getOperand(0)->getName());
                    } else {
#ifdef IDX_DEBUG
                        errs() << "Warning: PHI node without a name";
#endif
                    }
                } else if (isa<Instruction>(pInst->getOperand(0))) {
                    Instruction* tmp = dyn_cast<Instruction>(pInst->getOperand(0));
                    if (tmp->isBinaryOp()) {
                        expr.push_back(tmp->getOpcodeName());
                    } else if (isa<LoadInst>(tmp)) {
                        LoadInst* ld = dyn_cast<LoadInst>(tmp);
                        expr.push_back(ld->getOperand(0)->getName());
                    }
                }
            } else if (pInst->isBinaryOp()) {
                
                std::string op(pInst->getOpcodeName());
                //expr.push_back(op);
                
                if (pInst->getNumOperands() == 2) {
                    if (!isa<Instruction>(pInst->getOperand(0))) {
                        if(isa<Argument>(pInst->getOperand(0))) {
                            expr.push_back(pInst->getOperand(0)->getName().str());
                        } else if (isa<Constant>(pInst->getOperand(0))) {
                            Constant* cons = dyn_cast<Constant>(pInst->getOperand(0));
                            expr.push_back(cons->getUniqueInteger().toString(10, true));
                        } else {
#ifdef IDX_DEBUG
                            errs() << "Warning: Other operands\n";
#endif
                        }
                    } else if (isa<PHINode>(pInst->getOperand(0))) {
                        if (pInst->getOperand(0)->hasName()) {
                            expr.push_back(pInst->getOperand(0)->getName());
                        } else {
#ifdef IDX_DEBUG
                            errs() << "Warning: PHI node without a name";
#endif
                        }
                    } else if (isa<LoadInst>(pInst->getOperand(0))) {
                        LoadInst* ld = dyn_cast<LoadInst>(pInst->getOperand(0));
                        expr.push_back(ld->getOperand(0)->getName());
                    } else if (isa<Instruction>(pInst->getOperand(0))) {
                        Instruction* tmp = dyn_cast<Instruction>(pInst->getOperand(0));
                        if (tmp->isBinaryOp()) {
                            expr.push_back(tmp->getOpcodeName());
                        }
                    }
                    if (!isa<Instruction>(pInst->getOperand(1))) {
                        if(isa<Argument>(pInst->getOperand(1))) {
                            expr.push_back(pInst->getOperand(1)->getName().str());
                        } else if (isa<Constant>(pInst->getOperand(1))) {
                            Constant* cons = dyn_cast<Constant>(pInst->getOperand(1));
                            expr.push_back(cons->getUniqueInteger().toString(10, true));
                        } else {
#ifdef IDX_DEBUG
                            errs() << "Warning: Other operands\n";
#endif
                        }
                    } else if (isa<PHINode>(pInst->getOperand(1))) {
                        if (pInst->getOperand(1)->hasName()) {
                            expr.push_back(pInst->getOperand(1)->getName());
                        } else {
#ifdef IDX_DEBUG
                            errs() << "Warning: PHI node without a name";
#endif
                        }
                    } else if (isa<LoadInst>(pInst->getOperand(1))) {
                        LoadInst* ld = dyn_cast<LoadInst>(pInst->getOperand(1));
                        expr.push_back(ld->getOperand(0)->getName());
                    } else if (isa<Instruction>(pInst->getOperand(1))) {
                        Instruction* tmp = dyn_cast<Instruction>(pInst->getOperand(1));
                        if (tmp->isBinaryOp()) {
                            expr.push_back(tmp->getOpcodeName());
                        }
                    }
                    
                } else if (pInst->getNumOperands() == 1) {
                    if (!isa<Instruction>(pInst->getOperand(0))) {
                        if(isa<Argument>(pInst->getOperand(0))) {
                            expr.push_back(pInst->getOperand(0)->getName().str());
                        } else if (isa<Constant>(pInst->getOperand(0))) {
                            Constant* cons = dyn_cast<Constant>(pInst->getOperand(0));
                            expr.push_back(cons->getUniqueInteger().toString(10, true));
                        } else {
#ifdef IDX_DEBUG
                            errs() << "Warning: Other operands\n";
#endif
                        }
                    }   else if (isa<PHINode>(pInst->getOperand(0))) {
                        if (pInst->getOperand(0)->hasName()) {
                            expr.push_back(pInst->getOperand(0)->getName());
                        } else {
#ifdef IDX_DEBUG
                            errs() << "Warning: PHI node without a name";
#endif
                        }
                    } else if (isa<LoadInst>(pInst->getOperand(0))) {
                        LoadInst* ld = dyn_cast<LoadInst>(pInst->getOperand(0));
                        expr.push_back(ld->getOperand(0)->getName());
                    } else if (isa<Instruction>(pInst->getOperand(0))) {
                        Instruction* tmp = dyn_cast<Instruction>(pInst->getOperand(0));
                        if (tmp->isBinaryOp()) {
                            expr.push_back(tmp->getOpcodeName());
                        }
                    }
                }
            } else if (pInst->isShift()) {
#ifdef IDX_DEBUG
                errs() << "Error: shift instruction\n";
#endif
            } else if (isa<CallInst>(pInst)) {
                CallInst* callInst = dyn_cast<CallInst>(pInst);
                Value* valueCalled = callInst->getCalledValue();
                if (valueCalled == NULL) {
#ifdef IDX_DEBUG
                    errs() << "NULL value here\n";
#endif
                } else {
                    if (valueCalled->stripPointerCasts() != NULL) {
                        if (valueCalled->stripPointerCasts()->hasName()) {
                            
                            std::string dimension = "";
                            if (callInst->getNumArgOperands() == 1) {
                                if (isa<Constant>(callInst->getArgOperand(0))) {
                                    Constant* cons = dyn_cast<Constant>(callInst->getArgOperand(0));
                                    dimension = cons->getUniqueInteger().toString(10, true);
                                }
                            }
                            
                            if (valueCalled->stripPointerCasts()->getName() == "get_global_id") {
                                expr.push_back("gbid" + dimension);
                            }
                            if (valueCalled->stripPointerCasts()->getName() == "get_global_size") {
                                expr.push_back("gbs" + dimension);
                            }
                            if (valueCalled->stripPointerCasts()->getName() == "get_local_id") {
                                expr.push_back("lid" + dimension);
                            }
                            if (valueCalled->stripPointerCasts()->getName() == "get_local_size") {
                                expr.push_back("ls" + dimension);
                            }
                            if (valueCalled->stripPointerCasts()->getName() == "get_group_id") {
                                expr.push_back("gid" + dimension);
                            }
                        }
                    }
                }
            //} else if (isa<LoadInst>(pInst)) {
            //    expr.push_back(pInst->getOperand(0)->getName());
            //} else if (isa<AllocaInst>(pInst)) {
            //    expr.push_back(pInst->getName());
            } else {
#ifdef IDX_DEBUG
                errs() << "Warning: Other instructions\n";
#endif
            }
        }
        
        if (expr.size() == 0) {
#ifdef IDX_DEBUG
            return "No expression";
#endif
        }
        
#ifdef IDX_DEBUG
        for (unsigned long i = 0; i < expr.size(); i++) {
            errs() << "Prefix expr: " << expr[i] << "\n";
        }
        errs() << "\n";
#endif
        
        /* Prefix expression to infix expression */
        std::string expr_infix;
        
        while (expr.size() > 1) {
            long i;
#ifdef IDX_DEBUG
            for (i = 0; i < (long) expr.size(); i++) {
                errs() << "    " << expr[i];
            }
            errs() << "\n";
#endif
            for (i = expr.size() - 1; i >= 0; i--) {
                if (expr[i] == "add" || expr[i] == "sub" || expr[i] == "mul" || expr[i] == "div") {
                    break;
                }
            }
            
            if (expr[i] == "add") {
                expr[i] = "(" + expr[expr.size() - 2] + " + " + expr[expr.size() - 1] + ")";
            }
            if (expr[i] == "sub") {
                expr[i] = "(" + expr[expr.size() - 2] + " - " + expr[expr.size() - 1] + ")";
            }
            if (expr[i] == "mul") {
                expr[i] = "(" + expr[expr.size() - 2] + " * " + expr[expr.size() - 1] + ")";
            }
            if (expr[i] == "div") {
                expr[i] = "(" + expr[expr.size() - 2] + " / " + expr[expr.size() - 1] + ")";
            }
            expr.pop_back();
            expr.pop_back();
        }
        
        expr_infix = expr[0];
        
        return expr_infix;
    }
    
    /* Extracting array name from load/store instruction */
    std::string IndexAnalysis::getArrayName(GetElementPtrInst* inst) {
        std::string nameTmp = "";
        
        if (isa<LoadInst>(inst->getPointerOperand())) {
            /* Load instruction */
            LoadInst* ldTmp = dyn_cast<LoadInst>(inst->getPointerOperand());
            nameTmp = ldTmp->getOperand(0)->getName();
        } else {
            /* Store instruction */
            nameTmp = inst->getPointerOperand()->getName();
        }
        
        return nameTmp;
    }
    
    /* Scan all IR instructions in function to find load/stores accessing arrays */
    void IndexAnalysis::findAllArrayAccesses(Function &F) {
        
        /* Using instruction iterator to interating through all the IR instructions in the functions */
        for (inst_iterator it = inst_begin(F), eit = inst_end(F); it != eit; ++it) {
            
            if (isa<LoadInst>(*it)) {
                if (isa<GetElementPtrInst>(it->getOperand(0))) {
                    GetElementPtrInst* gepTmp = dyn_cast<GetElementPtrInst>(it->getOperand(0));
                    
                    /* Problem: how to figure out which operand is index for GEP? */
                    int OperandToAnalysis = 1;
                    
                    if (gepTmp->getNumOperands() < 2) {
#ifdef IDX_DEBUG
                        errs() << "ERROR: GEP does not have enough operands\n";
#endif
                    } else {
                        Instruction* accessInst = dyn_cast<Instruction>(&(*it));
                        std::string nameTmp = getArrayName(gepTmp);
                        
                        if (isa<Instruction>(gepTmp->getOperand(OperandToAnalysis))) {
                            Instruction* indexInst = dyn_cast<Instruction>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = computeExpression(indexInst);
                        } else if (isa<Argument>(gepTmp->getOperand(OperandToAnalysis))) {
                            Argument* argTmp = dyn_cast<Argument>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = argTmp->getName().str();
                        } else if (isa<Constant>(gepTmp->getOperand(OperandToAnalysis))) {
                            Constant* constTmp = dyn_cast<Constant>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = constTmp->getUniqueInteger().toString(10, true);
                        } else {
#ifdef IDX_DEBUG
                            errs() << "ERROR: other load address type\n";
#endif
                        }
                    }
                }
            }
            
            if (isa<StoreInst>(*it)) {
                if (isa<GetElementPtrInst>(it->getOperand(1))) {
                    GetElementPtrInst* gepTmp = dyn_cast<GetElementPtrInst>(it->getOperand(1));
                    
                    /* Problem: how to figure out which operand is index for GEP? */
                    int OperandToAnalysis = 1;
                    
                    if (gepTmp->getNumOperands() < 2) {
#ifdef IDX_DEBUG
                        errs() << "ERROR: GEP does not have enough operands\n";
#endif
                    } else {
                        Instruction* accessInst = dyn_cast<Instruction>(&(*it));
                        std::string nameTmp = getArrayName(gepTmp);
                        
                        if (isa<Instruction>(gepTmp->getOperand(OperandToAnalysis))) {
                            Instruction* indexInst = dyn_cast<Instruction>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = computeExpression(indexInst);
                        } else if (isa<Argument>(gepTmp->getOperand(OperandToAnalysis))) {
                            Argument* argTmp = dyn_cast<Argument>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = argTmp->getName().str();
                        } else if (isa<Constant>(gepTmp->getOperand(OperandToAnalysis))) {
                            Constant* constTmp = dyn_cast<Constant>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = constTmp->getUniqueInteger().toString(10, true);
                        } else {
#ifdef IDX_DEBUG
                            errs() << "ERROR: other store address type\n";
#endif
                        }
                    }
                }
            }
        }
        return;
    }
    
    /* Dump (array name, index expression) pairs for all references inside the function.
       Two hash tables are used to store array names and index expressions, the keys are load/store instruction pointers */
    void IndexAnalysis::dumpAllInfo() {
        errs() << "Array index info\n";
        for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            errs() << it->second << " " << arrayExpression[it->first] << "\n";
        }
        return;
    }
    
    /* Analysis main function: information are dumped as comments at the beginning of generated static sampling code */
    bool IndexAnalysis::runOnFunction(Function &F) {
        errs() << "\n /* Start to analysis array index\n";
        
        findAllArrayAccesses(F);
        
        dumpAllInfo();
        
        errs() << "\n Finish to analysis array index */ \n";
        
        return false;
    }
    
}

