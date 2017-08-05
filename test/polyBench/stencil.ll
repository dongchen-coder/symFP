; ModuleID = 'stencil.bc'
source_filename = "/Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/polyBench/stencil.c"
target datalayout = "e-m:o-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-apple-macosx10.12.0"

; Function Attrs: noinline nounwind ssp uwtable
define i32 @stencil(double* %a, double* %b) #0 {
entry:
  %a.addr = alloca double*, align 8
  %b.addr = alloca double*, align 8
  %i = alloca i32, align 4
  %j = alloca i32, align 4
  store double* %a, double** %a.addr, align 8
  store double* %b, double** %b.addr, align 8
  store i32 0, i32* %i, align 4
  br label %for.cond

for.cond:                                         ; preds = %for.inc31, %entry
  %0 = load i32, i32* %i, align 4
  %cmp = icmp slt i32 %0, 1025
  br i1 %cmp, label %for.body, label %for.end33

for.body:                                         ; preds = %for.cond
  store i32 0, i32* %j, align 4
  br label %for.cond1

for.cond1:                                        ; preds = %for.inc, %for.body
  %1 = load i32, i32* %j, align 4
  %cmp2 = icmp slt i32 %1, 1025
  br i1 %cmp2, label %for.body3, label %for.end

for.body3:                                        ; preds = %for.cond1
  %2 = load double*, double** %a.addr, align 8
  %3 = load i32, i32* %i, align 4
  %mul = mul nsw i32 1026, %3
  %4 = load i32, i32* %j, align 4
  %add = add nsw i32 %mul, %4
  %idxprom = sext i32 %add to i64
  %arrayidx = getelementptr inbounds double, double* %2, i64 %idxprom
  %5 = load double, double* %arrayidx, align 8
  %6 = load double*, double** %a.addr, align 8
  %7 = load i32, i32* %i, align 4
  %mul4 = mul nsw i32 1026, %7
  %8 = load i32, i32* %j, align 4
  %add5 = add nsw i32 %mul4, %8
  %add6 = add nsw i32 %add5, 1
  %idxprom7 = sext i32 %add6 to i64
  %arrayidx8 = getelementptr inbounds double, double* %6, i64 %idxprom7
  %9 = load double, double* %arrayidx8, align 8
  %add9 = fadd double %5, %9
  %10 = load double*, double** %a.addr, align 8
  %11 = load i32, i32* %i, align 4
  %mul10 = mul nsw i32 1026, %11
  %12 = load i32, i32* %j, align 4
  %add11 = add nsw i32 %mul10, %12
  %sub = sub nsw i32 %add11, 1
  %idxprom12 = sext i32 %sub to i64
  %arrayidx13 = getelementptr inbounds double, double* %10, i64 %idxprom12
  %13 = load double, double* %arrayidx13, align 8
  %add14 = fadd double %add9, %13
  %14 = load double*, double** %a.addr, align 8
  %15 = load i32, i32* %i, align 4
  %add15 = add nsw i32 %15, 1
  %mul16 = mul nsw i32 1026, %add15
  %16 = load i32, i32* %j, align 4
  %add17 = add nsw i32 %mul16, %16
  %idxprom18 = sext i32 %add17 to i64
  %arrayidx19 = getelementptr inbounds double, double* %14, i64 %idxprom18
  %17 = load double, double* %arrayidx19, align 8
  %add20 = fadd double %add14, %17
  %18 = load double*, double** %a.addr, align 8
  %19 = load i32, i32* %i, align 4
  %sub21 = sub nsw i32 %19, 1
  %mul22 = mul nsw i32 1026, %sub21
  %20 = load i32, i32* %j, align 4
  %add23 = add nsw i32 %mul22, %20
  %idxprom24 = sext i32 %add23 to i64
  %arrayidx25 = getelementptr inbounds double, double* %18, i64 %idxprom24
  %21 = load double, double* %arrayidx25, align 8
  %add26 = fadd double %add20, %21
  %22 = load double*, double** %b.addr, align 8
  %23 = load i32, i32* %i, align 4
  %mul27 = mul nsw i32 1026, %23
  %24 = load i32, i32* %j, align 4
  %add28 = add nsw i32 %mul27, %24
  %idxprom29 = sext i32 %add28 to i64
  %arrayidx30 = getelementptr inbounds double, double* %22, i64 %idxprom29
  store double %add26, double* %arrayidx30, align 8
  br label %for.inc

for.inc:                                          ; preds = %for.body3
  %25 = load i32, i32* %j, align 4
  %inc = add nsw i32 %25, 1
  store i32 %inc, i32* %j, align 4
  br label %for.cond1

for.end:                                          ; preds = %for.cond1
  br label %for.inc31

for.inc31:                                        ; preds = %for.end
  %26 = load i32, i32* %i, align 4
  %inc32 = add nsw i32 %26, 1
  store i32 %inc32, i32* %i, align 4
  br label %for.cond

for.end33:                                        ; preds = %for.cond
  ret i32 0
}

attributes #0 = { noinline nounwind ssp uwtable "correctly-rounded-divide-sqrt-fp-math"="false" "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-jump-tables"="false" "no-nans-fp-math"="false" "no-signed-zeros-fp-math"="false" "no-trapping-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "unsafe-fp-math"="false" "use-soft-float"="false" }

!llvm.module.flags = !{!0}
!llvm.ident = !{!1}

!0 = !{i32 1, !"PIC Level", i32 2}
!1 = !{!"clang version 4.0.0 (tags/RELEASE_400/final)"}
