//
//  idxAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/18/17.
//
//


#include "idxAnalysis.hpp"

using namespace llvm;

namespace idxAnalysis {
    
    char IndexAnalysis::ID = 0;
    static RegisterPass<IndexAnalysis> X("idxAnalysis", "Array index analysis Pass", false, false);
    
    IndexAnalysis::IndexAnalysis() : FunctionPass(ID), arrayName(), arrayExpression() {}
    
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
                errs() << "ERROR: Expression should not be that long\n";
                break;
            }
        }
        return instStack;
    }
    
    std::string IndexAnalysis::computeExpression(Instruction *inst) {
        
        inst->dump();
        
        std::vector<Instruction*> instStack = instStackInit(inst);
        
        if (instStack.size() == 0) {
            return "No expression";
        }
        
        /* iterate over the instructions stacked to get the prefix expression */
        
        std::vector<std::string> expr;
        for (std::vector<Instruction*>::iterator it = instStack.begin(), eit = instStack.end(); it != eit; ++it) {
            
            Instruction* pInst = *it;
            
            if (pInst->isCast()) {
                //stack[i]->getOperand(0)->dump();
                if (isa<PHINode>(pInst->getOperand(0))) {
                    if (pInst->getOperand(0)->hasName()) {
                        expr.push_back(pInst->getOperand(0)->getName());
                    } else {
                        errs() << "Warning: PHI node without a name";
                    }
                }
            } else if (pInst->isBinaryOp()) {
                
                std::string op(pInst->getOpcodeName());
                expr.push_back(op);
                
                if (pInst->getNumOperands() == 2) {
                    if (!isa<Instruction>(pInst->getOperand(0))) {
                        if(isa<Argument>(pInst->getOperand(0))) {
                            expr.push_back(pInst->getOperand(0)->getName().str());
                        } else if (isa<Constant>(pInst->getOperand(0))) {
                            Constant* cons = dyn_cast<Constant>(pInst->getOperand(0));
                            expr.push_back(cons->getUniqueInteger().toString(10, true));
                        } else {
                            errs() << "Warning: Other operands\n";
                        }
                    } else if (isa<PHINode>(pInst->getOperand(0))) {
                        if (pInst->getOperand(0)->hasName()) {
                            expr.push_back(pInst->getOperand(0)->getName());
                        } else {
                            errs() << "Warning: PHI node without a name";
                        }
                    }
                    if (!isa<Instruction>(pInst->getOperand(1))) {
                        if(isa<Argument>(pInst->getOperand(1))) {
                            expr.push_back(pInst->getOperand(1)->getName().str());
                        } else if (isa<Constant>(pInst->getOperand(1))) {
                            Constant* cons = dyn_cast<Constant>(pInst->getOperand(1));
                            expr.push_back(cons->getUniqueInteger().toString(10, true));
                        } else {
                            errs() << "Warning: Other operands\n";
                        }
                    } else if (isa<PHINode>(pInst->getOperand(1))) {
                        if (pInst->getOperand(1)->hasName()) {
                            expr.push_back(pInst->getOperand(1)->getName());
                        } else {
                            errs() << "Warning: PHI node without a name";
                        }
                    }
                    
                }
                
                if (pInst->getNumOperands() == 1) {
                    if (!isa<Instruction>(pInst->getOperand(0))) {
                        if(isa<Argument>(pInst->getOperand(0))) {
                            expr.push_back(pInst->getOperand(0)->getName().str());
                        } else if (isa<Constant>(pInst->getOperand(0))) {
                            Constant* cons = dyn_cast<Constant>(pInst->getOperand(0));
                            expr.push_back(cons->getUniqueInteger().toString(10, true));
                        } else {
                            errs() << "Warning: Other operands\n";
                        }
                    }   else if (isa<PHINode>(pInst->getOperand(0))) {
                        if (pInst->getOperand(0)->hasName()) {
                            expr.push_back(pInst->getOperand(0)->getName());
                        } else {
                            errs() << "Warning: PHI node without a name";
                        }
                    }
                }
            } else if (pInst->isShift()) {
                errs() << "Error: shift instruction\n";
            } else if (isa<CallInst>(pInst)) {
                CallInst* callInst = dyn_cast<CallInst>(pInst);
                Value* valueCalled = callInst->getCalledValue();
                if (valueCalled == NULL) {
                    errs() << "NULL value here\n";
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
            } else if (isa<AllocaInst>(pInst)) {
                expr.push_back(pInst->getName());
            } else {
                errs() << "Warning: Other instructions\n";
            }
        }
        
        if (expr.size() == 0) {
            return "No expression";
        }
        
        /* Prefix expression to infix expression */
        std::string expr_infix;
        
        int len = expr.size() - 1;
        for (int i = len; i >= 0; i--) {
            //errs() << expr[i] << " ";
            bool flag = false;
            if (expr[i] == "add") {
                expr[i] = "(" + expr[i+1] + " + " + expr[i+2] + ")";
                flag = true;
            }
            if (expr[i] == "sub") {
                expr[i] = "(" + expr[i+1] + " - " + expr[i+2] + ")";
                flag = true;
            }
            if (expr[i] == "mul") {
                expr[i] = "(" + expr[i+1] + " * " + expr[i+2] + ")";
                flag = true;
            }
            if (expr[i] == "div") {
                expr[i] = "(" + expr[i+1] + " / " + expr[i+2] + ")";
                flag = true;
            }
            
            if (flag == true) {
                for (int j = i+1; j <= len-2; j++) {
                    expr[j] = expr[j+2] ;
                }
                len = len - 2;
            }
        }
        
        expr_infix = expr[0];
        
        return expr_infix;
    }
    
    std::string IndexAnalysis::getArrayName(GetElementPtrInst* inst) {
        std::string nameTmp = "";
        
        if (isa<LoadInst>(inst->getPointerOperand())) {
            LoadInst* ldTmp = dyn_cast<LoadInst>(inst->getPointerOperand());
            nameTmp = ldTmp->getOperand(0)->getName();
        } else {
            nameTmp = inst->getPointerOperand()->getName();
        }
        
        return nameTmp;
    }
    
    void IndexAnalysis::findAllArrayAccesses(Function &F) {
        for (inst_iterator it = inst_begin(F), eit = inst_end(F); it != eit; ++it) {
            if (isa<LoadInst>(*it)) {
                if (isa<GetElementPtrInst>(it->getOperand(0))) {
                    GetElementPtrInst* gepTmp = dyn_cast<GetElementPtrInst>(it->getOperand(0));
                    
                    /* Problem: how to figure out which operand is index for GEP? */
                    int OperandToAnalysis = 1;
                    
                    if (gepTmp->getNumOperands() < 2) {
                        errs() << "ERROR: GEP does not have enough operands\n";
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
                            errs() << "ERROR: other load address type\n";
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
                        errs() << "ERROR: GEP does not have enough operands\n";
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
                            errs() << "ERROR: other store address type\n";
                        }
                    }
                }
            }
        }
        return;
    }
    
    void IndexAnalysis::dumpAllInfo() {
        errs() << "Array index info\n";
        for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            errs() << it->second << " " << arrayExpression[it->first] << "\n";
        }
        return;
    }
    
    bool IndexAnalysis::runOnFunction(Function &F) {
        errs() << "\nStart to analysis array index\n";
        
        F.dump();
        
        findAllArrayAccesses(F);
        
        dumpAllInfo();
        
        return false;
    }
    
}


